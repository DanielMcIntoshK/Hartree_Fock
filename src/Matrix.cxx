#include "Matrix.h"
#include <iostream>
#include <cmath>

std::vector<int> sortlist(std::vector<double> &sl){
	std::vector<int> sv;
	sv.resize(sl.size());

	for(int i = 0; i < sv.size(); i++){
		int place=0;
		for(int l = 0; l < i; l++){
			if(sl[l] <= sl[i]) place++;
		}
		for(int u=i+1;u<sv.size();u++){
			if(sl[u] < sl[i]) place++;

		}
		sv[place]=i;
	}

	std::vector<double> nl;
	nl.resize(sl.size());
	for(int i = 0; i < nl.size(); i++) nl[i]=sl[sv[i]];
	sl=nl;
	return sv;
}

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

void Matrix::printMatrix(){
	for(int r = 0; r < rows; r++){
		for(int c = 0; c < cols; c++){
			std::string instr="             ";
			std::string datastr=std::to_string(data[r*cols+c]);
			instr.replace(0,datastr.size()-1,datastr);
			std::cout << instr;
		}
		std::cout << std::endl;
	}
}

Matrix Matrix::matMul(Matrix& m1, Matrix& m2){
	if(m1.cols!=m2.rows){return Matrix(0,0);}
	Matrix nm(m1.rows, m2.cols);

	for(int r=0; r < nm.rows; r++){
		for(int c = 0; c < nm.cols; c++){
			double val=0.0;
			for(int i = 0; i < m1.cols; i++){
				val+= m1.data[r*m1.cols+i]*m2.data[i*m2.cols+c];
			}
			nm.data[r*nm.cols+c]=val;
		}
	}
	return nm;
}

Matrix Matrix::scaleMul(double s,Matrix & m){
	Matrix nm(m.rows,m.cols);
	for(int i = 0; i < m.data.size();i++){ nm.data[i]=m.data[i]*s;}
	return nm;
}

Matrix Matrix::transpose(){
	Matrix tm(cols,rows);
	for(int i = 0; i < rows;i++){
		for(int j=0;j<cols;j++){
			tm.setElement(j,i,data[i*cols+j]);
		}
	}
	return tm;
}

Matrix Matrix::diag(std::vector<double> eigenvals){
	int N=eigenvals.size();
	Matrix dm(N,N);
	for(int i = 0; i < N;i++){
		for(int j = 0; j < N;j++){
			if(i==j)dm.data[i*N+j]=std::pow(eigenvals[i],-.5);
			else dm.data[i*N+j]=0.0;
		}
	}
	return dm;
}

Matrix Matrix::matAdd(Matrix & m1, Matrix & m2){
	if(m1.rows!=m2.rows || m1.cols!=m2.cols) return Matrix(0,0);
	int N=m1.rows;
	int M=m1.cols;
	Matrix nm(N,M);
	for(int i = 0; i < N; i++){
		for(int j = 0; j < M; j++){
			int index = i*M+j;
			nm.data[index]=m1.data[index]+m2.data[index];
		}
	}
	return nm;
}

void Matrix::sort(std::vector<int> sortvals){
	std::vector<double> newdata;
	newdata.resize(data.size());

	for(int i = 0 ; i < rows; i++){
		for(int j = 0; j < cols;j++){
			newdata[i*cols+j]=data[i*cols+sortvals[j]];
		}
	}
	data=newdata;

}

void Matrix::identity(){
	for(int i = 0; i < rows; i++){
		for(int j = 0; j < cols; j++){
			if(i==j) data[i*cols+j]=1.0;
			else data[i*cols+j]=0.0;
		}
	}
}

EigenSolver::EigenSolver(Matrix & inMat, double _threshold):threshold{_threshold},oMat{inMat}{
	//oMat.printMatrix(8);
}

EigenSolver::EigenData EigenSolver::calculateEigens(){
	Matrix cMat=oMat;
	EigenData ed(oMat.rows);

	ed.eigenVecs.identity();
	for(int i = 0; i < ed.eigenVals.size();i++) ed.eigenVals[i]=0.0;

	int citer=0;
	int maxiter=cMat.rows*cMat.rows*cMat.rows;

	findMaxUpperTriangle(cMat);
	do{
		std::cout << citer << std::endl;
		jacobiRotate(cMat, ed);
		findMaxUpperTriangle(cMat);
		citer++;
	}while(std::fabs(cMat(maxr,maxc))>threshold && citer < maxiter);
	if(citer==maxiter) std::cout << "TIMEOUT\n";
	for(int n = 0; n < cMat.rows; n++) ed.eigenVals[n]=cMat(n,n);
	std::cout << "EIGENVALS: ";
	for(auto a: ed.eigenVals) std::cout << a << " ";
	std::cout << std::endl;
	ed.eigenVecs=ed.eigenVecs.transpose();
	return ed;
}

void EigenSolver::findMaxUpperTriangle(Matrix & cMat){
	maxr=0;
	maxc=0;
	double maxval=-1;
	for(int i = 0; i < cMat.rows-1; i++){
		for( int j = 1+i; j < cMat.cols; j++){
			if(std::fabs(cMat(i,j)) > maxval){
				maxval=std::fabs(cMat(i,j));
				maxr=i;
				maxc=j;
			}
		}
	}
}

void EigenSolver::jacobiRotate(Matrix & cMat, EigenSolver::EigenData & ed){
	double c, s;

	//std::cout << maxr << " " << maxc << " " << cMat(maxr,maxc) << std::endl;

	if(cMat(maxc,maxr) != 0.0){
		double tau, t;
		tau = (cMat(maxc,maxc)-cMat(maxr,maxr)) / (2.0 * cMat(maxc,maxr));
		if(tau > 0.0000001){
			t = 1.0/(tau+std::sqrt(1.0+tau*tau));
		}
		else{
			t = -1.0/(-tau+std::sqrt(1.0+tau*tau));
		}
		c=1.0/std::sqrt(1+t*t);
		s=c*t;
	}
	else{
		c=1.0;
		s=0.0;
	}


	std::vector<double> eigenVecsIRow,eigenVecsJRow;
	eigenVecsIRow.resize(cMat.rows);
	eigenVecsJRow.resize(cMat.rows);

	for(int j = 0; j < cMat.rows; j++){
		eigenVecsIRow[j]=ed.eigenVecs(maxr,j);
		eigenVecsJRow[j]=ed.eigenVecs(maxc,j);
		
		ed.eigenVecs.setElement(maxr,j,eigenVecsIRow[j]*c-eigenVecsJRow[j]*s);
		ed.eigenVecs.setElement(maxc,j,eigenVecsIRow[j]*s+eigenVecsJRow[j]*c);
	}

	double eigen_ii=cMat(maxr,maxr);
	double eigen_jj=cMat(maxc,maxc);
	double eigen_il, eigen_jk;

	//cMat.printMatrix();

	cMat.setElement(maxr,maxr,c*c*eigen_ii-2.0*s*c*cMat(maxr,maxc)+s*s*eigen_jj);
	cMat.setElement(maxc,maxc,s*s*eigen_ii+2.0*s*c*cMat(maxr,maxc)+c*c*eigen_jj);
	cMat.setElement(maxr,maxc,0.0);
	cMat.setElement(maxc,maxr,0.0);

	for(int l = 0; l < cMat.rows; l++){
		if(l != maxr && l != maxc){
			eigen_il=cMat(maxr,l);
			eigen_jk=cMat(maxc,l);

			//std::cout << eigen_il << " " << eigen_jk << std::endl;

			cMat.setElement(maxr,l,c*eigen_il-s*eigen_jk);
			cMat.setElement(l,maxr,cMat(maxr,l));

			cMat.setElement(maxc,l,s*eigen_il+c*eigen_jk);
			cMat.setElement(l,maxc,cMat(maxc,l));
		}
	}
	//cMat.printMatrix();
}
