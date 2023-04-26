#include "HF.h"
#include <iostream>
#include <cmath>

HF::HF(Matrix & _S, Matrix & _X, Matrix & _cH, list4D ld):
S{_S},X{_X},cH{_cH},eeList{ld}
{
	F=Matrix{S.rows,S.cols};
	P=Matrix{S.rows,S.cols};
	Xt=X.transpose();	
}

double HF::calculateEnergy(Molecule & mol, Matrix & P_init, int charge, bool verb){
	verbose=verb;

	//if(verbose){
	std::cout << "******************\n"<<
		"STARTING SFC CALC\n" <<
		"******************\n";
	//}
	electroncount=0;
	for(int i = 0; i < mol.atoms.size(); i++){
		electroncount += mol.atoms[i].atnum;
	}
	electroncount-=charge;
	
	dim=mol.orbitalcount;
	
	P=P_init;

	int stepCount = 1;

	do{
		std::cout << std::endl<<"SCF STEP " << stepCount++<<std::endl;
		HFStep();
		std::cout << "ENERGY DIFF: " << eDiff<<std::endl;
		std::cout << "TOTAL ENERGY: " << calculateEnergyTotal(mol)<<std::endl;
	}while(eDiff > 1.0E-12&&stepCount < 10000);

	if(stepCount==10000) std::cout << "\nSCF DID NOT CONVERGE\n";
	else std::cout << "\nSCF CONVERGED\n";
	return calculateEnergyTotal(mol);
}

void HF::HFStep(){
	if(verbose){
	std::cout << "\nDENSITY MATRIX:\n";
	P.printMatrix();
	std::cout << std::endl;
	}

	if(verbose) std::cout << "Making FOCK Matrix\n";
	computeInteraction(P);
	if(verbose) std::cout << "CALC ee interaction\n";

	F = Matrix::matAdd(cH,eeI);
	
	if(verbose){
	std::cout << "FOCK MATRIX\n";
	F.printMatrix();

	std::cout << "Transform Fock Matrix\n";
	}

	Matrix Fx=Matrix::matMul(Xt,F);
	Fx=Matrix::matMul(Fx,X);

	if(verbose){
	std::cout << std::endl << "FOCK MATRIX TRANSFORMED\n";
	Fx.printMatrix();
	std::cout << std::endl;
	

	std::cout << "Calculating EigenVectors\n";
	}

	EigenSolver es(Fx,0.000001);
	EigenSolver::EigenData ed = es.calculateEigens();
	
	if(verbose){
	std::cout << "\nEIGENVALS: ";
	for(auto &a: ed.eigenVals) std::cout << a << " ";
	std::cout << std::endl;

	std::cout << "\nEIGENVECS:\n";
	ed.eigenVecs.printMatrix();
	}

	std::vector<int> sl=sortlist(ed.eigenVals);
	ed.eigenVecs.sort(sl);
	
	if(verbose){
	std::cout << "\nEIGENVALS: ";
	for(auto &a: ed.eigenVals) std::cout << a << " ";
	std::cout << std::endl;

	std::cout << "\nEIGENVECS:\n";
	ed.eigenVecs.printMatrix();
	}

	Matrix C = Matrix::matMul(X,ed.eigenVecs);

	if(verbose){
	std::cout << "\nCOEFFICIENT MATRIX\n";
	C.printMatrix();
	}

	Matrix P_New=computeNewDensity(C);

	computeDensityDifference(P,P_New);

	P=P_New;
}

void HF::computeInteraction(Matrix & P_i){
	eeI=Matrix(dim,dim);

	for(int i=0; i < dim; i++){
		for(int j = 0; j < dim; j++){
			double sub=0.0;
			for(int k = 0; k < dim; k++){
				for(int l = 0; l < dim; l++){
					sub+= P_i(k,l)*(eeList[i][j][k][l]-0.5*eeList[i][l][k][j]);	
				}
			}
			eeI.setElement(i,j,sub);
		}
	}
}

Matrix HF::computeNewDensity(Matrix & C){
	Matrix newDen(C.rows, C.cols);

	for( int i = 0; i < C.rows; i++){
		for( int j = 0; j < C.cols; j++){
			for(int k = 0; k < (int)(electroncount/2);k++){
				newDen.data[i*newDen.cols+j]+=2*C(i,k)*C(j,k);
			}
		}
	}
	return newDen;
}

void HF::computeDensityDifference(Matrix & P_Old, Matrix & P_New){
	double delta = 0.0;

	int dim = P_Old.rows;

	for(int i = 0; i < dim; i++){
		for(int j =0; j < dim; j++){
			double diff = P_Old(i,j)-P_New(i,j);
			delta+=diff*diff;
		}
	}
	eDiff= std::sqrt(delta/4.0);
}

double HF::calculateEnergyTotal(Molecule & mol){
	return calculateENuc(mol)+calculateEElectron();
}

double HF::calculateENuc(Molecule & mol){
	double En=0.0;
	for(int i = 0; i < mol.atoms.size(); i++){
		for(int j = i+1; j < mol.atoms.size(); j++){
			Atom & a1=mol.atoms[i], & a2 = mol.atoms[j];

			double dprod=0.0;
			for(int i = 0; i < 3; i++){
				double diff=a1.pos[i]-a2.pos[i];
				dprod+=diff*diff;
			}
			double Z1=a1.atnum, Z2 = a2.atnum;
			En+=Z1*Z2/(std::sqrt(dprod));
		}
	}	
	return En;
}

double HF::calculateEElectron(){
	double E = 0.0;
	for(int i = 0; i < P.rows;i++){
		for(int j=0;j<P.cols;j++){
			E+=0.5*P(i,j)*(cH(i,j)+F(i,j));
		}
	}
	return E;
}


