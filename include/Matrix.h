#ifndef MATRIX__H_
#define MATRIX__H_
#include <vector>

class Matrix{
public:
	Matrix(int r, int c);

	void setData(std::vector<double> d);
	void setElement(int i, int c, double d);

	//double operator[](int r, int c){return data[r*cols+c];}
	double operator()(int r, int c){return data[r*cols+c];}

	static Matrix matMul(Matrix & m1, Matrix &m2);

	std::vector<double> data;
	int rows, cols;
};


#endif

