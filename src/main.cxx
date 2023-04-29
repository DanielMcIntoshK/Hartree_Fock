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
	mol.loadMolecule("Molecules/ethanol.xyz");
	//Molecule H2O;
	//mol.atoms=std::vector<Atom>{Atom(1,{1.4305507,0,0}),Atom(1,{-1.4305507,0,0}),Atom(8,{0,1.1072514,0})};
	//mol.init("basissets/def2-svp.1.orca");
	
	mol.init("basissets/sto-3g.orca");
	if(procno==0){
		std::cout << "MOLECULE: " << mol.orbitalcount << std::endl;
		mol.printMol();
	}

	double ti=MPI_Wtime();

	bool verbose=false;

	DistributedMatrix S=computeSMatrixPar(mol);
	if(verbose){
		if(S.procno==0)std::cout << "OVERLAP MATRIX\n";
		S.printMatrix();
	}
	double ot=MPI_Wtime();

	DistributedEigenSolver des(S, 0.000000001);
	DistributedEigenSolver::EigenData ded = des.calculateEigens();
	
	double eigt=MPI_Wtime();

	if(verbose){
	if(S.procno==0){
		std::cout << "EIGEN VALS: ";
		for(auto a: ded.eigenVals) std::cout << a << " ";
		std::cout << std::endl;
		}
		//std::cout << "EIGENVECS:\n";
		ded.eigenVecs.printMatrix();
	}
	for(int i = 0; i < ded.eigenVals.size(); i++){
		ded.eigenVals[i]=1/std::sqrt(ded.eigenVals[i]);
	}
	DistributedMatrix diag=DistributedMatrix::diag(ded.eigenVals);

	DistributedMatrix X=DistributedMatrix::matMul(ded.eigenVecs,diag);
	DistributedMatrix Xt=X.transpose();
	
	if(verbose){
		if(X.procno==0)std::cout << "\nTRANSFORM MATRIX:\n";
		X.printMatrix();
		if(X.procno==0)std::cout << std::endl;
	}
	double prehamt=MPI_Wtime();

	DistributedMatrix cH=computeCoreHamiltonianMatrixPar(mol);

	double posthamt=MPI_Wtime();

	if(verbose){
	if(cH.procno==0)std::cout << "CORE HAMILTONIAN\n";
	cH.printMatrix();
	if(cH.procno==0)std::cout << std::endl;
	}

	DistributedMatrix ld=computeEEMatriciesPar(mol);
	std::vector<double> ldMat=ld.gatherMat();

	double eecomputet=MPI_Wtime();

	DistributedMatrix P_init(S.rows,S.cols,MPI_COMM_WORLD);

	HFPar hartreefockSolver(S,X,cH,ld,ldMat);
	hartreefockSolver.calculateEnergy(mol,P_init,0,verbose);

	double end=MPI_Wtime();

	if(procno==0){
		std::cout << "Program exited after " << end-ti << " seconds.\n";
		std::cout <<"*******\n"<<"OVERALAP MATRIX: " << ot-ti << std::endl <<
			"Eigenvec algorithm: " << eigt-ot << std::endl <<
			"Core Hamiltonian: " << posthamt-prehamt<<std::endl <<
			"EE integral list: " << eecomputet-posthamt << std::endl <<
			"SCF Converge: " << end-eecomputet << std::endl;
	}

	MPI_Finalize();

	return 0;

}
