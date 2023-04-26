#ifndef MATRIX__H_
#define MATRIX__H_
#include <vector>

std::vector<int> sortlist(std::vector<double> &sl);

class Matrix{
public:
	Matrix():rows{0},cols{0}{}
	Matrix(int r, int c);

	void setData(std::vector<double> d);
	void setElement(int r, int c, double d);

	//double operator[](int r, int c){return data[r*cols+c];}
	
	double operator()(int r, int c){return data[r*cols+c];}

	void printMatrix();

	static Matrix matMul(Matrix & m1, Matrix &m2);
	static Matrix scaleMul(double s, Matrix &m);
	static Matrix diag(std::vector<double> eigenvals);
	static Matrix matAdd(Matrix & m1, Matrix &m2);

	void sort(std::vector<int> sortvals);

	Matrix transpose();
	void identity();

	std::vector<double> data;
	int rows, cols;
};

class EigenSolver{
public:
	struct EigenData{
		EigenData(int s):eigenVecs{s,s}{eigenVals.resize(s);}
		Matrix eigenVecs;
		std::vector<double> eigenVals;
	};

	EigenSolver(Matrix & inMat, double _threshold);
	
	EigenData calculateEigens();

	double threshold;
	Matrix oMat;
private:
	int maxr, maxc;

	void findMaxUpperTriangle(Matrix & cMat);
	void jacobiRotate(Matrix & cmat, EigenData & ed);
};


#endif

