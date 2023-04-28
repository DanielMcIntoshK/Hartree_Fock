#include "Atom.h"
#include "Matrix.h"
#include "Integrals.h"
#include "MyMath.h"
#include "HF.h"
#include "DistributedMatrix.h"
#include <iostream>
#include <string>
#include <cmath>
#include <mpi.h>

int main(int argc, char ** argv){
	MPI_Init(&argc, &argv);
	MPI_Comm wcomm = MPI_COMM_WORLD;

	int nprocs, procno;
	MPI_Comm_size(wcomm, &nprocs);
	MPI_Comm_rank(wcomm, &procno);

	Molecule mol;
	mol.loadMolecule("Molecules/waterJank.xyz");
	//Molecule H2O;
	//mol.atoms=std::vector<Atom>{Atom(1,{1.4305507,0,0}),Atom(1,{-1.4305507,0,0}),Atom(8,{0,1.1072514,0})};
	//mol.init("basissets/def2-svp.1.orca");
	
	mol.init("basissets/sto-3g.orca");
	std::cout << "MOLECULE: " << mol.orbitalcount << std::endl;
	mol.printMol();

	bool verbose=true;

	DistributedMatrix S=computeSMatrixPar(mol);
	Matrix Sm=computeSMatrix(mol);
	if(verbose){
		if(S.procno==0)std::cout << "OVERLAP MATRIX\n";
		S.printMatrix();
	}

	DistributedEigenSolver des(S, 0.000000001);
	DistributedEigenSolver::EigenData ded = des.calculateEigens();
	
	EigenSolver es(Sm, 0.000001);
	EigenSolver::EigenData ed = es.calculateEigens();

	if(verbose){
	if(S.procno==0){
		std::cout << "EIGEN VALS: ";
		for(auto a: ded.eigenVals) std::cout << a << " ";
		std::cout << std::endl;
		}
		ded.eigenVecs.printMatrix();
		DistributedMatrix t=ded.eigenVecs.transpose();
		t.printMatrix();
	}
	for(int i = 0; i < ded.eigenVals.size(); i++){
		ded.eigenVals[i]=1/std::sqrt(ded.eigenVals[i]);
	}
	DistributedMatrix diag=DistributedMatrix::diag(ded.eigenVals);
	diag.printMatrix();

	Matrix mdiag=Matrix::diag(ed.eigenVals);
	Matrix Xm = Matrix::matMul(ed.eigenVecs,mdiag);
	
	
	DistributedMatrix X=DistributedMatrix::matMul(ded.eigenVecs,diag);
	//DistributedMatrix eigenVecstrans=ded.eigenVecs.transpose();
	//X=DistributedMatrix::matMul(X,eigenVecstrans);
	DistributedMatrix Xt=X.transpose();
	Xt.printMatrix();

	
	if(verbose){
		if(X.procno==0)std::cout << "\nTRANSFORM MATRIX:\n";
		X.printMatrix();
		if(X.procno==0)std::cout << std::endl;
	}


	DistributedMatrix test = DistributedMatrix::matMul(Xt,S);
	test=DistributedMatrix::matMul(test,X);
	if(test.procno==0)std::cout << "TEST SHOULD BE IDENT\n";
	test.printMatrix();

	
	DistributedMatrix cH=computeCoreHamiltonianMatrixPar(mol);
	Matrix cHm=computeCoreHamiltonianMatrix(mol);

	if(verbose){
	if(cH.procno==0)std::cout << "CORE HAMILTONIAN\n";
	cH.printMatrix();
	if(cH.procno==0)std::cout << std::endl;
	}

	DistributedMatrix ld=computeEEMatriciesPar(mol);
	list4D ldm = computeEEMatricies(mol);

	std::vector<double> ldmat=ld.gatherMat();
	if(S.procno==0){
	for(int p = 0; p < ldmat.size();p++){
		int dim = S.rows;
		int i=(p/(dim*dim*dim))%dim,
		    j=(p/(dim*dim))%dim,
		    k=(p/(dim))%dim,
		    l=p%dim;

		std::cout << "("<<i+1<<","<<j+1<<","<<k+1<<","<<l+1<<") "<<ldmat[p]<<" "<<ldm[i][j][k][l]<<std::endl;
	}
	}

	DistributedMatrix P_init(S.rows,S.cols,MPI_COMM_WORLD);
	Matrix P_initm(Sm.rows, Sm.cols);

	if(S.procno==0){
		HF hartreefockSolver(Sm,Xm,cHm,ldm);
		hartreefockSolver.calculateEnergy(mol,P_initm,0,verbose);
	}

	HFPar hartreefockSolver(S,X,cH,ld);
	hartreefockSolver.calculateEnergy(mol,P_init,0,verbose);

	MPI_Finalize();

	return 0;

}
