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

	bool verbose=false;

	Matrix S=computeSMatrix(mol);
	
	DistributedMatrix dm(S.rows, S.cols, wcomm);

	for(int i = 0; i < dm.data.size(); i++){
		int index = i+dm.startpos;
		int gr = dm.i2r(index);
		int gc = dm.i2c(index);

		dm.data[i]=S(gr,gc);
	}

	dm.printMatrix();

	EigenSolver no(S,0.000001);
	if(dm.procno==0){
		std::cout <<"SERIALEIGEN\n"<<std::endl;
		EigenSolver::EigenData nod=no.calculateEigens();
		std::cout << "DONE\n";
	}
	MPI_Barrier(wcomm);

	DistributedEigenSolver des(dm, 0.000001);
	DistributedEigenSolver::EigenData ded = des.calculateEigens();

	if(dm.procno==0){
		std::cout << "EIGEN VALS: ";
		for(auto a: ded.eigenVals) std::cout << a << " ";
		std::cout << std::endl;
	}
	ded.eigenVecs.printMatrix();
	
	MPI_Barrier(wcomm);

	MPI_Finalize();

	return 0;



	if(verbose){
	std::cout << "S MATRIX:\n";
	S.printMatrix();
	}

	Matrix cH=computeCoreHamiltonianMatrix(mol);

	if(verbose){
	std::cout << "CORE HAMILTONIAN\n";
	cH.printMatrix();
	}

	list4D ld=computeEEMatricies(mol);
	

	for(int i = 0; i < ld.size(); i++){
	for(int j = 0; j < ld[i].size(); j++){
	for(int k = 0; k < ld[i][j].size(); k++){
	for(int l = 0; l < ld[i][j][k].size();l++){
		if(std::isnan(ld[i][j][k][l])){
			std::cout << "NAN: " << i << " " << j << " " << k << " " << l << std::endl;
		}	
	}}}}

	EigenSolver es(S, 0.000001);
	EigenSolver::EigenData ed = es.calculateEigens();

	/*
	std::vector<int> sl{0,3,1,4,2,5,6};
	ed.eigenVecs.sort(sl);
	std::vector<double> sv;
	sv.resize(sl.size());
	for(int i = 0;i < sv.size();i++) sv[i]=ed.eigenVals[sl[i]];
	ed.eigenVals=sv;

	
	sl=std::vector<int>{1,2,3};
	for(auto a: sl){
		for(int r=0;r<ed.eigenVecs.cols;r++){
			std::cout << r << " " << a << " " << ed.eigenVecs(r,a)<<std::endl;
			ed.eigenVecs.setElement(r,a,-ed.eigenVecs(r,a));
		}
	}
	
	std::cout << "\nEIGENVALS: ";
	for(auto a: ed.eigenVals) std::cout << a << " ";
	std::cout << std::endl;
	std::cout << "EIGENVECS:\n";
	ed.eigenVecs.printMatrix();
	std::cout << std::endl;
	*/
	Matrix diag=Matrix::diag(ed.eigenVals);
	//Matrix eigenVecinv=ed.eigenVecs.transpose();

	//Matrix temp=Matrix::matMul(diag,eigenVecinv);
	Matrix X = Matrix::matMul(ed.eigenVecs,diag);
	if(verbose){
	std::cout << std::endl;
	std::cout << "TRANSFORM MATRIX\n";
	X.printMatrix();
	}

	Matrix Xt=X.transpose();
	
	if(verbose){
	Matrix com = Matrix::matMul(Xt,S);
	com=Matrix::matMul(com,X);
	std::cout << "COMBINED TEST\n";
	com.printMatrix();
	}

	Matrix P_init(S.rows,S.cols);

	HF hartreefockSolver(S,X,cH,ld);
	hartreefockSolver.calculateEnergy(mol,P_init,0,verbose);

	return 0;
}
