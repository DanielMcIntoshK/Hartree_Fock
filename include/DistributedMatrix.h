#ifndef DISTRIBUTEDMATRIX__H_
#define DISTRIBUTEDMATRIX__H_
#include <mpi.h>
#include <vector>


class DistributedMatrix{
public:
	DistributedMatrix():DistributedMatrix{MPI_COMM_WORLD}{}
	DistributedMatrix(MPI_Comm cm);
	DistributedMatrix(int r, int c, MPI_Comm cm);

	int rc2i(int r, int c){return r*cols+c;}
	int i2r(int i){return i/cols;}
	int i2c(int i){return i%cols;}
	int i2proc(int i);

	static DistributedMatrix matMul(DistributedMatrix & m1, DistributedMatrix & m2);
	static DistributedMatrix diag(std::vector<double> eigenvals);
	static DistributedMatrix matAdd(DistributedMatrix & m1, DistributedMatrix & m2);

	void sort(std::vector<int> sortvals);

	std::vector<double> gatherMat();

	DistributedMatrix transpose();
	void identity();

	void printMatrix();

	MPI_Comm comm;

	std::vector<double> data;
	int rows, cols, size, blocksize;

	int startpos;

	int blocksizebase, overflow;

	int procno, nprocs;
private:
};

class DistributedEigenSolver{
public:
	struct EigenData{
		EigenData(int s, MPI_Comm cm):eigenVecs{s,s,cm}{eigenVals.resize(s);}
		DistributedMatrix eigenVecs;
		std::vector<double> eigenVals;
	};
	struct valloc{
		double val;
		int loc;
	};

	DistributedEigenSolver(DistributedMatrix & inMat, double _t);

	EigenData calculateEigens();

	double threshold;
	DistributedMatrix oMat;
private:
	MPI_Comm comm;

	valloc maxvl;
	double maxvalTrans;
	int maxr, maxc;

	double diagr, diagc;

	void findMaxUpperTriangle(DistributedMatrix &cMat);
	void jacobiRotate(DistributedMatrix & cmat, EigenData & ed);
};

/*
class DistributedList4D{
public:
	DistributedList4D(MPI_Comm cm);

	std::vector<double> data;

	
private:

};
*/
#endif

