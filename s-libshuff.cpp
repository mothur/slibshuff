
//
// S-LIBSHUFF: A program for pairwise comparison of sequence libraries.
// Copyright (C) 2003, 2004, 2005, 2006 Bret Larget and Patrick D. Schloss 
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//

#define VERSION "1.22"
#define DATE "April 16, 2006"

// Version history:
// Version 1.22
//   outputs delCxy values
// Version 1.21
//   expanded input file formats to include mega3 formatted distance matrices
// Version 1.2
// 	 corrected a small bug that has gone largely undetected that behaves poorly when numGroups>9 or so
// 	 added output file feature that contains coverage values for input distance matrix
// 	 released initial version of manual
// 	 Released November 5, 2004
// Version 1.1
// 	 added ability to use lower triangular matrix for input file
// 	 Released May 25, 2004
// Version 1.0
//   fixed error in sCalculate (forgetting j!=i)
//   changed the output format
//   Released February 27, 2004
// Version 0.3 updated by Bret Larget
//   fixed bug?
//   Released January 16, 2004
// Version 0.2 written by Bret Larget and Patrick D. Schloss
//   Released January 9, 2003
// Version 0.1 written by Bret Larget and Patrick D. Schloss
//   Released December 28, 2003

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <list>
#include <cmath>

using namespace std;

// *******************************************************
//
// Random number generator modified from R, version 1.8.1
// http://www.r-project.org/
// R-1.8.1/src/nmath/standalone/sunif.c
//

static unsigned int I1=1234, I2=5678;

double runif(void)
{
    I1= 36969*(I1 & 0177777) + (I1>>16);
    I2= 18000*(I2 & 0177777) + (I2>>16);
    return ((I1 << 16)^(I2 & 0177777)) * 2.328306437080797e-10; /* in [0,1) */
}

void printSeeds(const unsigned int I1,const unsigned int I2)
{
  cout << "Seeds: " << I1 << " " << I2 << endl;
}
  
// *******************************************************

void swap(int& i,int& j)
{
  int t = i;  i = j;  j = t;
}

void swap(vector<int>& x,int a,int b)
{
  int t = x[a]; x[a] = x[b]; x[b] = t;
}

void swap(vector<double>& x,int a,int b)
{
  double t = x[a]; x[a] = x[b]; x[b] = t;
}

void sort(vector<double> const& x,vector<int>& y,int a,int b)
{
  if(a>=b)
    return;
  swap(y,a,(a+b)/2);
  int c = a;
  for(int i=a+1; i<=b; i++)
    if(x[y[i]] < x[y[a]])
      swap(y,++c,i);
  swap(y,a,c);
  sort(x,y,a,c-1);
  sort(x,y,c+1,b);
}

void sort(vector<double>& x,int a,int b)
{
  if(a>=b)
    return;
  swap(x,a,(a+b)/2);
  int c = a;
  for(int i=a+1; i<=b; i++)
    if(x[i] < x[a])
      swap(x,++c,i);
  swap(x,a,c);
  sort(x,a,c-1);
  sort(x,c+1,b);
}

void sort(vector<double>& x)
{
  sort(x,0,x.size()-1);
}

class DistanceMatrix {
public:
  DistanceMatrix(){};
  double calculate(int,int);
  void calculateAll();
  void calculateFamilyErrorRate(void);
  void calcNX(int,vector<int>&);
  void calcX(int,vector<double>&);
  void calcNXY(int,int,vector<int>&);
  void calcXY(int,int,vector<double>&);
  void dCalculateAll();
  double dCalculate(int,int);
  double findAllPValues(void);
  double findMinimumPValue(void);
  double getDeltaXY(int i,int j,int p){ return deltaXY[i][j][p]; }
  int getFamily(){ return family; }
  int getN(){ return n; }
  int getNumFamilyPerms(){ return numFamilyPerms; }
  int getNumGroups(){ return numGroups; }
  int getNumPerms(){ return numPerms; }
  double getSavedDeltaXY(int i,int j){ return savedDeltaXY[i][j]; }
  int getSummationSize(){ return summationSize; };
  void initialize(int,char*[]);
  void initializeDeltaXY(void);
  void initializeGroups();
  void initializePValue(void);
  void minimumX(int, int, vector<vector<double> >&);
  void minimumXY(int, int, vector<vector<vector<double> > >&);
  void print(const vector<int>);
  void printDeltaXY(int);
  void printGroupSizes();
  void printGroups();
  void printPValues();
  void printdelCxyValues();
  void randomizeAllGroups();
  void randomizeGroups(int,int);
  void read(int);
  void resetGroup(int);
  void resetPValueCounts(void);
  void saveAllGroups();
  void saveGroup(int);
  void setCutoffs(int);
  void setDeltaXY(int i,int j,int p,double dxy){ deltaXY[i][j][p] = dxy; }
  void setGroupSizes();
  void printcoverageD(void);
  void printcoverageS(void);
  void coverageS(vector<double>, vector<double>);
  double sCalculate(int,int);
  void sCalculateAll();
private:
	void read_phylip(istream&, int);
	void read_mega(istream&);
  int n; // total number of sequences
  int numGroups; // number of libraries
  int matrix; // 1=square matrix, 2=lower triangular matrix
  int numPerms; // number of permutations (default is 10000)
  int numFamilyPerms; // number of permutations (default is 10000)
  int discrete; // 0: use the integral form of the calculation; 1: use the discrete form of the calculation
  int family; // 0: do not compute family-wise error rate; 1: do compute family-wise error rate
  int summationSize; // number of distance categories used in discrete calculation
  int printout;
  vector<vector<double> > d; // n by n pairwise distance matrix
  vector<int> sizes; // vector of numGroup group sizes
  vector<vector<int> > groups; // groups[i][j] is the index of the jth sequence in group i
  vector<vector<int> > savedGroups; // savedGroups[i][j] is the index of the jth sequence in savedGroup i
  vector<double> cutoffs; // endpoints for discrete calculation
  vector<vector<vector<double> > > deltaXY; // deltaXY[x][y][p] is the test statistic for test of group x versus y for permutation p
  vector<vector<double> > savedDeltaXY; // test statistic for each test for orginal groups
  vector<vector<double> > coverage;
  int coverageindex;
  vector<vector<int> > pValueCounts; // count of extreme test statistics for each test
  string fileName;
};

/**************************************************************************************************/

void get_comment(istream& f, char begin, char end)
{
	char d=f.get();
	while(d != end){	d = f.get();	}
	d = f.peek();
}	

/**************************************************************************************************/

void DistanceMatrix::read(int matrix)
{
	ifstream f(fileName.c_str());
	if(!f) {
		cerr << "Error: Could not open " << fileName << endl;
		exit(1);
	}
	char test = f.peek();
	
	if(test == '#'){
		read_mega(f);
	}
	else{
		read_phylip(f, matrix);
	}
}

/**************************************************************************************************/

void DistanceMatrix::read_mega(istream& f)
{
	int count1 = 0;
	int count2 = 0;
	n = 0;
	get_comment(f, '#', '\n');

	char test = f.peek();

	while(test == '!'){								//get header comments
		get_comment(f, '!', ';');
		while(isspace(test=f.get()))		{;}
		f.putback(test);
		test = f.peek();
	}
	while(test != '\n'){							//get sequence names
		get_comment(f, '[', ']');
		char d = f.get();
		d = f.get();
		if(d == '#'){
			string name;
			f >> name;
			n++;
//			name_list.push_back(name);
			while(isspace(test=f.get()))		{;}
			f.putback(test);
		}
		else{
			break;
		}
	}
	d.resize(n);
	for(int i=0;i<n;i++){		d[i].resize(n);		}

	d[0][0] = 0.0000;
	count2++;
	get_comment(f, '[', ']');
	count1++;
	for(int i=1;i<n;i++){
		get_comment(f, '[', ']');
		d[i][i]=0.0000;
		count2++;
		for(int j=0;j<i;j++){
			f >> d[i][j];
			count2++;
			if (d[i][j] == -0.0000)
				d[i][j] = 0.0000;
			d[j][i]=d[i][j];
			count2++;
		}
		count1++;
	}
	cout << "Read a total of " << count2 << " distances, " << count1 << " individuals." << endl;
}

/**************************************************************************************************/

void DistanceMatrix::read_phylip(istream& f, int matrix)
{
	int count1=0,count2=0;
	
	f >> n;
	d.resize(n);
	
	if(matrix==1){
		for(int i=0;i<n;i++)
			d[i].resize(n);
		for(int i=0;i<n;i++){
			string name;
			f >> name;
			for(int j=0;j<n;j++) {	
				f >> d[i][j];
				count2++;    
			}
			count1++;
		}
	}
	else if(matrix==2){
		for(int i=0;i<n;i++){
			d[i].resize(n);
		}
		string name;
		d[0][0] = 0.0000;
		f >> name;
		count1++;
		for(int i=1;i<n;i++){
			f >> name;
			d[i][i]=0.0000;
			for(int j=0;j<i;j++){
				f >> d[i][j];
				count2++;
				if (d[i][j] == -0.0000)
					d[i][j] = 0.0000;
				d[j][i]=d[i][j];
				count2++;
			}
			count1++;
		}
		count2+=count1;
	}
	cout << "Read a total of " << count2 << " distances, " << count1 << " individuals." << endl;
}

void DistanceMatrix::print(const vector<int> rows)
{
  int m = rows.size();
  for(int i=0;i<m;i++) {
    for(int j=0;j<m;j++)
      cout << setw(9) << d[rows[i]][rows[j]];
    cout << endl;
  }
}

void DistanceMatrix::setGroupSizes()
{
  cout << "Enter number of groups: ";
  cin >> numGroups;
  cout << endl;
  if( numGroups<1 || numGroups > n) {
//  if( numGroups<1 || numGroups >> n) {
    cerr << "Error: Number of Groups must be in the range [1.." << n << "]." << endl;
    exit(1);
  }
  sizes.resize(numGroups);

  int sum = 0;

  cout << "Enter the group sizes summing to " << n << " (separated by spaces). ";
  for(int i=0;i<numGroups;i++) {
    cin >> sizes[i];
    sum += sizes[i];
    if( sizes[i]<1 ) {
      cerr << "Error: Group sizes must be positive." << endl;
      exit(1);
    }
  }
  if(sum != n) {
    cerr << "Error: Group sizes total is " << sum << " and should be " << n << "." << endl;
    exit(1);
  }
}

void DistanceMatrix::printGroupSizes()
{
  for(int i=0;i<numGroups;i++)
    cout << setw(5) << sizes[i];
  cout << endl;
}

void DistanceMatrix::initializeGroups()
{
  groups.resize(numGroups);
  savedGroups.resize(numGroups);
  for(int i=0;i<numGroups;i++) {
    groups[i].resize(sizes[i]);
    savedGroups[i].resize(sizes[i]);
  }
  int index=0;
  for(int i=0;i<numGroups;i++)
    for(int j=0;j<sizes[i];j++)
      savedGroups[i][j] = groups[i][j] = index++;
}


void DistanceMatrix::randomizeGroups(int i,int j)
{
  int nv = sizes[i]+sizes[j];
  vector<int> v(nv);
  int index=0;
  for(int k=0;k<sizes[i];k++)
    v[index++] = groups[i][k];
  for(int k=0;k<sizes[j];k++)
    v[index++] = groups[j][k];
  for(int k=nv-1;k>0;k--) {
    int z = (int)((k+1)*runif());
    swap(v[z],v[k]);
  }
  index=0;
  for(int k=0;k<sizes[i];k++)
    groups[i][k]=v[index++];
  for(int k=0;k<sizes[j];k++)
    groups[j][k]=v[index++];
}

void DistanceMatrix::randomizeAllGroups()
{
  vector<int> p(n);
  int index=0;
  for(int i=0;i<numGroups;i++)
    for(int j=0;j<sizes[i];j++)
      p[index++] = groups[i][j];

  for(int j=n-1;j>0;j--) {
    int i = (int)((j+1)*runif());
    swap(p[i],p[j]);
  }

  index=0;
  for(int i=0;i<numGroups;i++)
    for(int j=0;j<sizes[i];j++)
      groups[i][j] = p[index++];
}


void DistanceMatrix::resetGroup(int i)
{
  for(int k=0;k<sizes[i];k++)
    groups[i][k] = savedGroups[i][k];
}

void DistanceMatrix::saveGroup(int i)
{
  for(int k=0;k<sizes[i];k++)
    savedGroups[i][k] = groups[i][k];
}

void DistanceMatrix::saveAllGroups(void)
{
  for(int i=0;i<numGroups;i++)
    saveGroup(i);
}

void DistanceMatrix::setCutoffs(int s)
{
  cutoffs.resize(s);
  double c = 0.5/(s-1);
  cutoffs[0] = 0.0;
  for(int k=1;k<s;k++){
    cutoffs[k] = cutoffs[k-1] + c;
  }
}

void DistanceMatrix::initializeDeltaXY()
{
  deltaXY.resize(numGroups);
  savedDeltaXY.resize(numGroups);
  for(int i=0;i<numGroups;i++) {
    deltaXY[i].resize(numGroups);
    savedDeltaXY[i].resize(numGroups);
    for(int j=0;j<numGroups;j++)
      deltaXY[i][j].resize(numPerms);
  }
}

void DistanceMatrix::initializePValue()
{
  pValueCounts.resize(numGroups);
  for(int i=0;i<numGroups;i++) {
    pValueCounts[i].resize(numGroups);
    for(int j=0;j<numGroups;j++)
      pValueCounts[i][j] = 0;
  }
}

void DistanceMatrix::resetPValueCounts()
{
  for(int i=0;i<numGroups;i++)
    for(int j=0;j<numGroups;j++)
      pValueCounts[i][j] = 0;
}

void DistanceMatrix::calcNX(int x,vector<int>& nx)
{
  int m = sizes[x];
  for(int k=0;k<summationSize;k++)
    nx[k] = 0;

  for(int i=0;i<m;i++) {
    double minX = 1.0;
    for(int j=0;j<m;j++)
      if(j != i) {
	double dx = d[groups[x][i]][groups[x][j]];
	if(dx < minX)
	  minX = dx;
      }
    int k=0;
    while( (cutoffs[k] < minX) && (k < summationSize) )
      nx[k++]++;
  }  
}

void DistanceMatrix::calcX(int x, vector<double>& cx)
{
  vector<int> nx(summationSize);
  calcNX(x,nx);
  for(int k=0;k<summationSize;k++)
    cx[k] = 1.0 - (double)nx[k]/(double)sizes[x];
}

void DistanceMatrix::calcNXY(int x,int y,vector<int>& nxy)
{
  for(int k=0;k<summationSize;k++)
    nxy[k] = 0;

  for(int i=0;i<sizes[x];i++) {
    double minX = 1.0;
    for(int j=0;j<sizes[y];j++) {
      double dxy = d[groups[x][i]][groups[y][j]];
      if(dxy < minX)
	minX = dxy;
    }
    int k=0;
    while( (cutoffs[k] < minX) && (k < summationSize) )
      nxy[k++]++;
  }
}

void DistanceMatrix::calcXY(int x,int y,vector<double>& cxy)
{
  vector<int> nxy(summationSize);
  calcNXY(x,y,nxy);
  for(int k=0;k<summationSize;k++)
    cxy[k] = 1.0 - (double)nxy[k]/(double)sizes[x];
}

double DistanceMatrix::dCalculate(int i,int j)
{
  if(i==j)
    return 0.0;
  vector<double> cx(summationSize);
  vector<double> cxy(summationSize);
  calcX(i,cx);
  calcXY(i,j,cxy);
 	double x = 0;
 
	if(printout == 1){
		coverage[coverageindex].resize(summationSize);
		coverage[coverageindex+1].resize(summationSize);
		
		for(int i=0;i<summationSize;i++){
  			coverage[coverageindex][i] = cx[i];
  			coverage[coverageindex+1][i] = cxy[i];
		}
		coverageindex += 2;
	}
  
  double sum = 0.0;
  for(int k=0;k<summationSize;k++) {
    double diff = cx[k] - cxy[k];
    sum += diff*diff;
  }
  return sum;
}

void DistanceMatrix::printcoverageD(void)
{
	string outputfile = fileName + ".coverage";

	ofstream f(outputfile.c_str(), ios::trunc);
	f.setf(ios::fixed, ios::floatfield);
	f.setf(ios::showpoint);

	f << setw(8) << "Dist";

	for(int i=0;i<numGroups;i++){
		f << setw(7) << "C" << i;
		for(int j=0;j<numGroups;j++){
			if(i!=j){
				f << setw(6) << "C" << i << j;
			}
		}
	}
	f << endl;
	for(int i=0;i<summationSize;i++){						//iterate across the distances
		f << setw(8) << setprecision(4) << coverage[0][i];	//print the distances

		int index = 1;
		while(index < coverage.size()){
			f << setw(8) << setprecision(4) << coverage[index][i];
			f << setw(8) << setprecision(4) << coverage[index+1][i];
			int index2;	
			for(index2=3;index2<2*(numGroups-1);index2+=2){
				f << setw(8) << setprecision(4) << coverage[index+index2][i];
			}
			index+=index2-1;
		}
		f << endl;
	}
	f.close();
}

void DistanceMatrix::printcoverageS(void)
{
	string outputfile = fileName + ".coverage";

	ofstream f(outputfile.c_str(), ios::trunc);
	f.setf(ios::fixed, ios::floatfield);
	f.setf(ios::showpoint);

	int maxlength = 0;
	for(int i=0;i<coverage.size();i++){
		if(coverage[i].size() > maxlength){
			maxlength = coverage[i].size();
		}
	}

	for(int i=0;i<numGroups;i++){
		for(int j=0;j<numGroups;j++){
			if(i!=j){
				f << setw(8) << "Dist" << setw(7) << "C" << i << setw(6) << "C" << i << j;
			}
		}
	}
	f << endl;

	for(int j=0;j<maxlength;j++){
		for(int i=0;i<coverage.size();i++){
			if(j+1>coverage[i].size()){
				f << setw(8) << "";
				continue;
			}
			else{
				f << setw(8) << setprecision(4) << coverage[i][j];
			}
		}
		f << endl;
	}
}

double DistanceMatrix::sCalculate(int x,int y)
{
  vector<double> minX(sizes[x]);
  vector<double> minXY(sizes[x]);
  for(int i=0;i<sizes[x];i++) {
    minX[i] = ( sizes[x] > 1 ? ( i==0 ? d[groups[x][0]][groups[x][1]] : d[groups[x][i]][groups[x][0]]) : 0.0);
    for(int j=0;j<sizes[x];j++) {
      if(j!=i) {
	double dx = d[groups[x][i]][groups[x][j]];
	if(dx < minX[i])
	  minX[i] = dx;
      }
    }
    minXY[i] = d[groups[x][i]][groups[y][0]];
    for(int j=0;j<sizes[y];j++) {
      double dxy = d[groups[x][i]][groups[y][j]];
      if(dxy < minXY[i])
	minXY[i] = dxy;
    }
  }
  sort(minX);
  sort(minXY);

  if(printout == 1){
	  coverageS(minX, minXY);
  }
  
  double sum = 0.0,t=0.0;
  int ix=0,iy=0;
  while( (ix < sizes[x]) && (iy < sizes[x]) ) {
    double h = ix-iy;
    if(minX[ix] < minXY[iy]) {
      sum += (minX[ix] - t)*h*h;
      t = minX[ix++];
    }
    else {
      sum += (minXY[iy] - t)*h*h;
      t = minXY[iy++];
    }
  }
  if(ix < sizes[x]) {
    while(ix < sizes[x]) {
      double h = ix-iy;
      sum += (minX[ix] - t)*h*h;
      t = minX[ix++];
    }
  }
  else {
    while(iy < sizes[x]) {
      double h = ix-iy;
      sum += (minXY[iy] - t)*h*h;
      t = minXY[iy++];
    }
  }
  return sum;// / (double)sizes[x];
}

void DistanceMatrix::coverageS(vector<double> X, vector<double> XY)
{
	int x_size = X.size();
	int xy_size = XY.size();
	
	{	//Get list of unique distances
		vector<double> tempdists;
		for(int i=0;i<x_size;i++){
			tempdists.push_back(X[i]);
		}
		for(int i=0;i<xy_size;i++){
			tempdists.push_back(XY[i]);
		}
		sort(tempdists.begin(),tempdists.end());
	
		coverage[coverageindex].push_back(tempdists[0]);
		for(int i=1;i<tempdists.size();i++){
			if(tempdists[i] != tempdists[i-1]){
				coverage[coverageindex].push_back(tempdists[i]);
			}
		}
	}
	double xbump = 1.0/(double)x_size;
	double xybump = 1.0/(double)xy_size;

	map<double,double> x_coverage;
	x_coverage[X[0]] = xbump;	
	for(int i=1;i<x_size;i++){
		x_coverage[X[i]] = xbump+x_coverage[X[i-1]];
	}

	map<double,double> xy_coverage;
	xy_coverage[XY[0]] = xybump;
	for(int i=1;i<xy_size;i++){
		xy_coverage[XY[i]] = xybump+xy_coverage[XY[i-1]];
	}

	int distances = coverage[coverageindex].size();
	coverage[coverageindex+1].resize(distances);
	coverage[coverageindex+2].resize(distances);

	for(int i=0;i<distances;i++){
		if(x_coverage[coverage[coverageindex][i]]){
			coverage[coverageindex+1][i] = x_coverage[coverage[coverageindex][i]];
		}
		else{
			coverage[coverageindex+1][i] = coverage[coverageindex+1][i-1];
		}
	
		if(xy_coverage[coverage[coverageindex][i]]){
			coverage[coverageindex+2][i] = xy_coverage[coverage[coverageindex][i]];
		}
		else{
			coverage[coverageindex+2][i] = coverage[coverageindex+2][i-1];
		}
	}
	coverageindex+=3;
}


double DistanceMatrix::calculate(int i,int j)
{
  if(discrete)
     return dCalculate(i,j);
  else
    return sCalculate(i,j);
}

void DistanceMatrix::dCalculateAll()
{
  for(int i=0;i<numGroups;i++)
    for(int j=0;j<numGroups;j++)
      if(i!=j)
		savedDeltaXY[i][j] = dCalculate(i,j);
      else
		savedDeltaXY[i][j] = 0.0;
}

void DistanceMatrix::sCalculateAll()
{
  for(int i=0;i<numGroups;i++)
    for(int j=0;j<numGroups;j++)
      if(i!=j)
		savedDeltaXY[i][j] = sCalculate(i,j);
      else
		savedDeltaXY[i][j] = 0.0;
}

void DistanceMatrix::calculateAll()
{
  printout = 1;
  if(discrete){
	coverageindex = 1;
	coverage.resize(2*numGroups*(numGroups-1)+1);
	coverage[0] = cutoffs;
	dCalculateAll();
	printcoverageD();
  }
  else{
	coverageindex = 0;
    coverage.resize(3*numGroups*(numGroups-1));
	sCalculateAll();
	printcoverageS();
  }
  printout = 0;
}

void DistanceMatrix::printDeltaXY(int p)
{
  cout << endl << "DeltaXY" << endl << endl;
  for(int i=0;i<numGroups;i++) {
    for(int j=0;j<numGroups;j++)
      cout << setprecision(4) << setw(7) << deltaXY[i][j][p];
    cout << endl;
  }
  cout << endl;
}

void DistanceMatrix::printPValues()
{
  cout << endl;
  cout << "   |      Y" << endl;
  cout << " X |";
  for(int j=0;j<numGroups;j++)
    cout << setw(7) << j+1;
  cout << endl;
  cout << "----";
  for(int j=0;j<numGroups;j++)
    cout << "-------";
  cout << endl;

  for(int i=0;i<numGroups;i++) {
    cout << setw(2) << i+1 << " |";
    for(int j=0;j<numGroups;j++)
      cout << setprecision(4) << setw(7) << (double)pValueCounts[i][j]/(double)numPerms;
    cout << endl;
  }
  cout << "----";
  for(int j=0;j<numGroups;j++)
    cout << "-------";
  cout << endl << endl;
}

void DistanceMatrix::printdelCxyValues()
{
  cout << endl;
  cout << "   |      Y" << endl;
  cout << " X |";
  for(int j=0;j<numGroups;j++)
    cout << setw(7) << j+1;
  cout << endl;
  cout << "----";
  for(int j=0;j<numGroups;j++)
    cout << "-------";
  cout << endl;

  for(int i=0;i<numGroups;i++) {
    cout << setw(2) << i+1 << " |";
    for(int j=0;j<numGroups;j++){
		if(i==j){	cout << setprecision(1) << setw(7) << 0.0000;						}
		else	{	cout << setprecision(1) << setw(7) << (double)savedDeltaXY[i][j];	}
	}
    cout << endl;
  }
  cout << "----";
  for(int j=0;j<numGroups;j++)
    cout << "-------";
  cout << endl << endl;
}


double DistanceMatrix::findAllPValues(void)
{
  int incr = numPerms/10;
  double smallestP=1.0;
  cout << endl << "// Beginning Calculation." << endl << endl;
  cout << "                   % Complete" << endl;
  cout << " Groups | 10 20 30 40 50 60 70 80 90 100" << endl;
  cout << "----------------------------------------" << endl;
  for(int i=0;i<numGroups-1;i++) {
    for(int j=i+1;j<numGroups;j++) {
      int check = incr;
      cout << setw(3) << i+1 << setw(4) << j+1 << " |";
      cout.flush();
      for(int p=0;p<numPerms;p++) {
	randomizeGroups(i,j);
	double dxy = calculate(i,j);
	setDeltaXY(i,j,p,dxy);
	if(deltaXY[i][j][p] >= savedDeltaXY[i][j])
	  pValueCounts[i][j]++;
	dxy = calculate(j,i);
	setDeltaXY(j,i,p,dxy);
	if(deltaXY[j][i][p] >= savedDeltaXY[j][i])
	  pValueCounts[j][i]++;
	if(p==check) {
	  cout << "***";
	  cout.flush();
	  check += incr;
	}
      }
      double pValue = (double)pValueCounts[i][j] / (double)numPerms;
      if(pValue < smallestP)
	smallestP = pValue;
      pValue = (double)pValueCounts[j][i] / (double)numPerms;
      if(pValue < smallestP)
	smallestP = pValue;
      resetGroup(i);
      resetGroup(j);
      cout << "****" << endl;
    }
  }
  cout << "----------------------------------------" << endl;
  cout << endl << "// dCxy Values" << endl;
  printdelCxyValues();

  cout << endl << "// P-values" << endl;
  printPValues();
  cout << "----------------------------------------" << endl;
  
  cout.flush();
  return smallestP;
}

double DistanceMatrix::findMinimumPValue()
{
  int minCount=numPerms;
  for(int i=0;i<numGroups-1;i++) {
    for(int j=i+1;j<numGroups;j++) {
      for(int p=0;p<numPerms;p++) {
	randomizeGroups(i,j);
	double dxy = calculate(i,j);
	setDeltaXY(i,j,p,dxy);
	if(deltaXY[i][j][p] >= savedDeltaXY[i][j])
	  pValueCounts[i][j]++;
	dxy = calculate(j,i);
	setDeltaXY(j,i,p,dxy);
	if(deltaXY[j][i][p] >= savedDeltaXY[j][i])
	  pValueCounts[j][i]++;
	if(pValueCounts[i][j]>=minCount && pValueCounts[j][i]>=minCount)
	  break;
      }
      if(pValueCounts[i][j] < minCount)
	minCount = pValueCounts[i][j];
      if(pValueCounts[j][i] < minCount)
	minCount = pValueCounts[j][i];
      resetGroup(i);
      resetGroup(j);
    }
  }
  return (double)minCount / (double)numPerms;
}

void usageError(char *name)
{
  cerr << "Usage: " << name << " [-s1 seedOne] [-s2 seedTwo] [-n numPermutations] [-d] [-z size] [-e numFamilyPermutations] <file>"
       << endl;
  cerr << " Options:" << endl;
  cerr << " -s1: first seed for random number generator, any positive integer (default is 1234)." << endl;
  cerr << " -s2: second seed for random number generator, any positive integer (default is 5678)." << endl;
  cerr << " -n: number of permutations for each permutation test (default is 10000)." << endl;
  cerr << " -d: use the discrete approximation for calculation (default is integral form)." << endl;
  cerr << " -z: size of increment for discrete calculation (only used if -d set; default is 0.01)." << endl;
  cerr << " -e: number of complete permutations to estimate family-wise error rate (default is 0)." << endl;
  cerr << " -l: lower triangular distance matrix (default is a square matrix)." << endl;
  cerr << " <file> is a symmetric distance matrix as produced by PHYLIP." << endl;
  exit(1);
}

void DistanceMatrix::initialize(int argc, char *argv[])
{
  cout.setf(ios::fixed, ios::floatfield);
  cout.setf(ios::showpoint);
  cerr.setf(ios::fixed, ios::floatfield);
  cerr.setf(ios::showpoint);

  cout << "// S-LIBSHUFF version " << setprecision(1) << VERSION << " (" << DATE << ")" << endl;
  cout << "// Written by Bret Larget and Pat Schloss." << endl << endl;

  cout << "// Call:";
  for(int i=0; i<argc; i++)
    cout << " " << argv[i];
  cout << endl;
  cout.flush();

  char **p;
  double z = 0.01;
  discrete = 0;
  family = 0;
  matrix = 1;
  numPerms = 10000;
  numFamilyPerms = 0;

  if(argc > 1) {
	  for(p=argv+1;p<argv+argc;p++) {
		  if(strcmp(*p,"-s1")==0) {
			  if(++p>=argv+argc)
				  usageError(argv[0]);
			  istringstream f(*p);
			  if(!(f >> I1))
				  usageError(argv[0]);
			  if(I1==0) {
				  cerr << "Error: random seeds must be positive." << endl;
				  usageError(argv[0]);
			  }
		  }
		  else if(strcmp(*p,"-s2")==0) {
			  if(++p>=argv+argc)
				  usageError(argv[0]);
			  istringstream f(*p);
			  if(!(f >> I2))
				  usageError(argv[0]);
			  if(I2==0) {
				  cerr << "Error: random seeds must be positive." << endl;
				  usageError(argv[0]);
			  }
		  }
		  else if(strcmp(*p,"-n")==0) {
			  if(++p>=argv+argc)
				  usageError(argv[0]);
			  istringstream f(*p);
			  if(!(f >> numPerms))
				  usageError(argv[0]);
		  }
		  else if(strcmp(*p,"-d")==0) {
			  if(p>=argv+argc)
				  usageError(argv[0]);
			  discrete=1;
		  }
		  else if(strcmp(*p,"-z")==0) {
			  if(++p>=argv+argc)
				  usageError(argv[0]);
			  istringstream f(*p);
			  if(!(f >> z))
				  usageError(argv[0]);
			  if(z < 1.0e-6) {
				  cerr << "Error: Size (-z) must be positive." << endl;
				  exit(1);
			  }
		  }
		  else if(strcmp(*p,"-e")==0) {
			  if(++p>=argv+argc)
				  usageError(argv[0]);
			  family=1;
			  istringstream f(*p);
			  if(!(f >> numFamilyPerms))
				  usageError(argv[0]);
		  }
		  else if(strcmp(*p,"-l")==0) {
			  if(p>=argv+argc)
				  usageError(argv[0]);
			  matrix=2;
		  }
		  else {
			  istringstream f(*p);
			  if(!(f >> fileName))
				  usageError(argv[0]);
		  }
	  }
	  read(matrix);
  }
  else
	  usageError(argv[0]);

  setGroupSizes();
  initializeGroups();
  initializeDeltaXY();
  initializePValue();
  summationSize = (int)(0.50 / z)+2;
  setCutoffs(summationSize);
  calculateAll();
}

void DistanceMatrix::calculateFamilyErrorRate()
{
  cout << endl << "// Minimum P-values for Family-wise Error Rate" << endl << endl;
  vector<double> minPValue(numFamilyPerms);
  int incr = numFamilyPerms/10;
  cout << "  Beginning Calculation." << endl << endl;
  cout << "         % Complete" << endl;
  cout << " 10 20 30 40 50 60 70 80 90 100" << endl;
  cout << "-------------------------------" << endl;
  int check = incr;
  for(int p=0;p<numFamilyPerms;p++) {
    randomizeAllGroups();
    saveAllGroups();
    calculateAll();
    resetPValueCounts();
    minPValue[p] = findMinimumPValue();
    if(p==check) {
      cout << "***";
      cout.flush();
      check += incr;
    }
  }
  cout << "***" << endl << "-------------------------------" << endl << endl;
  sort(minPValue);
  //  for(int p=0;p<numFamilyPerms;p++)
  //    cout << setprecision(4) << minPValue[p] << endl;
  double mean=0;
  for(int p=0;p<numFamilyPerms;p++)
    mean += minPValue[p];
  mean /= numFamilyPerms;
  double sd=0;
  for(int p=0;p<numFamilyPerms;p++) {
    double x=minPValue[p] - mean;
    sd += x*x;
  }
  sd = sqrt(sd / numFamilyPerms);
  double me = 2 * sd / sqrt((double)numFamilyPerms);
  cout << endl << "// 95% confidence interval for mean minimum P-Value is approximately:" << endl;
  cout << "    " << mean - me << " to " << mean + me;
  cout << "    " << "(Theory predicts " << 1.0 / (1.0 + numGroups*(numGroups-1)) << ".)" << endl;
  int count1=0;
  int count5=0;
  // Assumes sorted p-values
  while(minPValue[count5]<0.01) {
    count1++;
    count5++;
  }
  while(minPValue[count5]<0.05)
    count5++;
  cout << "    " << count1 << " of " << numFamilyPerms << " randomizations have a minimal p-value less than 0.01." << endl;
  cout << "    " << count5 << " of " << numFamilyPerms << " randomizations have a minimal p-value less than 0.05." << endl;
}

int main(int argc, char *argv[])
{
  DistanceMatrix d;
  d.initialize(argc,argv);
  double smallestP = d.findAllPValues();
  int m = d.getNumGroups()*(d.getNumGroups()-1);
  double a = 1.0 - pow(0.95,1.0/m);
  double b = 1.0 - pow(0.99,1.0/m);
  double c = 1.0 - pow(1.0-smallestP,m);
  cout << "// Correction for Multiple Comparisons (" << m << " tests)" << endl << endl;
  cout << "Family-wise error rate | Minimum p-value" << endl;
  cout << "----------------------------------------" << endl;
  if(c > 0.05) {
    cout << "   [ " << setprecision(4) << setw(6) << c << " ]          |";
    cout << "   [ " << setprecision(4) << setw(6) << smallestP << " ]" << endl;
  }
  cout << setprecision(4) << setw(11) << 0.05 << "            |";
  cout << setprecision(4) << setw(11) << a << endl;
  if(c > 0.01 && c < 0.05) {
    cout << "   [ " << setprecision(4) << setw(6) << c << " ]          |";
    cout << "   [ " << setprecision(4) << setw(6) << smallestP << " ]" << endl;
  }
  cout << setprecision(4) << setw(11) << 0.01 << "            |";
  cout << setprecision(4) << setw(11) << b << endl;
  if(c < 0.01) {
    cout << "   [ " << setprecision(4) << setw(6) << c << " ]          |";
    cout << "   [ " << setprecision(4) << setw(6) << smallestP << " ]" << endl;
  }
  int numPerms = d.getNumPerms();
  cout << endl << "// Monte Carlo Error ( " << numPerms << " permutations)" << endl << endl;
  cout << "95% margin of error for minimum p-value is "
       << 2*sqrt(smallestP * (1.0-smallestP)/(double)(numPerms));
  cout << endl;

  if(d.getFamily())
    d.calculateFamilyErrorRate();
  cout << endl;
  return 0;
}
