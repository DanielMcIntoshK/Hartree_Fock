#include "Atom.h"
#include "Matrix.h"
#include "Integrals.h"
#include "MyMath.h"
#include "HF.h"
#include <iostream>
#include <string>

int main(int argc, char ** argv){
	Molecule H2O;
	H2O.atoms=std::vector<Atom>{Atom(1,{1.4305507,0,0}),Atom(1,{-1.4305507,0,0}),Atom(8,{0,1.1072514})};
	H2O.init("basissets/sto-3g.orca");
	std::cout << "H2O: " << H2O.orbitalcount << std::endl;
	
	Matrix S=computeSMatrix(H2O);
	
	std::cout << "S MATRIX:\n";
	S.printMatrix();

	Matrix cH=computeCoreHamiltonianMatrix(H2O);

	std::cout << "CORE HAMILTONIAN\n";
	cH.printMatrix();

	list4D ld=computeEEMatricies(H2O);
	
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
	std::cout << std::endl;
	std::cout << "TRANSFORM MATRIX\n";
	X.printMatrix();

	Matrix Xt=X.transpose();
	Matrix com = Matrix::matMul(Xt,S);
	com=Matrix::matMul(com,X);
	std::cout << "COMBINED TEST\n";
	com.printMatrix();

	Matrix P_init(S.rows,S.cols);

	HF hartreefockSolver(S,X,cH,ld);
	hartreefockSolver.calculateEnergy(H2O,P_init,0,false);

	return 0;
}
