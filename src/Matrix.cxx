#include "Matrix.h"

Matrix::Matrix(int r, int c):rows{r},cols{c}{
	data.resize(rows*cols);
	for(auto& a:data){a=0;}
}

void Matrix::setData(std::vector<double> d){
	if(d.size() == rows*cols) data=d;
}

void Matrix::setElement(int r, int c, double d){
	data[r*cols+c]=d;
}

//Matrix& Matrix::operator*(Matrix& m){
	
//}

Matrix Matrix::matMul(Matrix& m1, Matrix& m2){
	if(m1.cols!=m2.rows){return Matrix(0,0);}
	Matrix nm(m1.rows, m2.cols);

	for(int r=0; r < nm.rows; r++){
		for(int c = 0; c < nm.cols; c++){
			double val=0.0;
			for(int i = 0; i < nm.cols; i++){
				val+= m1.data[r*m1.cols+i]*m2.data[i*m2.cols+c];
			}
			nm.data[r*nm.cols+c]=val;
		}
	}
	return nm;
}
