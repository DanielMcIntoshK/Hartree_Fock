#include "HF.h"
#include <iostream>
#include <cmath>

HFPar::HFPar(DistributedMatrix & _S, DistributedMatrix & _X, DistributedMatrix & _cH,
		DistributedMatrix & ld):S{_S},X{_X},cH{_cH},eeList{ld}{
	F=DistributedMatrix(S.rows,S.cols,S.comm);
	P=DistributedMatrix(S.rows,S.cols,S.comm);
	Xt=X.transpose();
}

double HFPar::calculateEnergy(Molecule & mol, DistributedMatrix & P_init, int charge, bool verb){
	verbose=verb;

	if(S.procno==0){
		std::cout << 
			"*****************\n" <<
			"STARTING SFC CALC\n" <<
			"*****************\n";
	}
	electroncount=0;
	for(int i = 0; i < mol.atoms.size(); i++){
		electroncount += mol.atoms[i].atnum;
	}
	electroncount-=charge;

	if(S.procno==0)std::cout << "ELECTRONS: " << electroncount<<std::endl;

	dim = mol.orbitalcount;

	P=P_init;

	int stepCount = 1;

	do{
		if(S.procno==0) std::cout << std::endl << "SCF STEP " << stepCount++ << std::endl;
		HFStep();
		double totalE = calculateEnergyTotal(mol);
		if(S.procno==0) std::cout << "ENERGY DIFF: " << eDiff << std::endl << "TOTAL ENERGY: " << totalE << std::endl;
	}while(eDiff > 1.0E-12 && stepCount < 10000);
	if(stepCount>=10000) std::cout << "\nSCF DID NOT CONVERGE\n";
	else std::cout << "\nSCF CONVERGED\n";

	return calculateEnergyTotal(mol);
}

void HFPar::HFStep(){
	
	if(S.procno==0){
		std::cout << "COMPUTING INTERACTION\n";
	}
	P.printMatrix();

	computeInteraction(P);
	
	if(S.procno==0){
		std::cout << "COMPUTING FOCK MATRIX\n";
	}
	F=DistributedMatrix::matAdd(cH,eeInteract);
	if(S.procno==0)std::cout << "EEINTERACT\n";
	eeInteract.printMatrix();
	F.printMatrix();
	if(S.procno==0){
		std::cout << "Transforming FOCK To ORTHOGANAL BASIS SET\n";
	}

	DistributedMatrix Fx = DistributedMatrix::matMul(Xt,F);
	Fx=DistributedMatrix::matMul(Fx,X);
	if(S.procno==0){
		std::cout << "Solving EigenVecs of Fock\n";
	}

	Fx.printMatrix();

	DistributedEigenSolver es(Fx,0.000001);
	DistributedEigenSolver::EigenData ed = es.calculateEigens();
	if(S.procno==0){
		std::cout << "EGENVALS: ";
		for(auto a: ed.eigenVals) std::cout << a << " ";
		std::cout << std::endl<<"EIGENVECS: \n";
	}
	ed.eigenVecs.printMatrix();
	if(S.procno==0){
		std::cout << "Sorting EigenVecs\n";
	}

	
	std::vector<int> sl=sortlist(ed.eigenVals);
	ed.eigenVecs.sort(sl);
	if(S.procno==0){
		std::cout << "Calculating C\n SORTED\n";
	}
	ed.eigenVecs.printMatrix();
	
	
	DistributedMatrix C = DistributedMatrix::matMul(X,ed.eigenVecs);
	if(S.procno==0){
		std::cout << "Computing New Density\n";
	}
	C.printMatrix();
	
	DistributedMatrix P_New=computeNewDensity(C);
	if(S.procno==0){
		std::cout << "Computing Density Difference\n";
	}
	P_New.printMatrix();

	computeDensityDifference(P,P_New);

	P=P_New;
	
}

void HFPar::computeInteraction(DistributedMatrix & Pc){
	eeInteract=DistributedMatrix(dim, dim,Pc.comm);
	int procno=Pc.procno;
	for(int i = 0; i < dim; i++){
		for(int j = 0; j < dim; j++){
			double sub=0.0;
			int dataproc=eeInteract.i2proc(eeInteract.rc2i(i,j));
			for(int k = 0; k < dim; k++){
				for(int l = 0; l < dim;l++){
					int indexij=i*dim*dim*dim+j*dim*dim+k*dim+l;
					int indexil=i*dim*dim*dim+l*dim*dim+k*dim+l;
					
					int indexkl=Pc.rc2i(k,l);

					int ijproc=eeList.i2proc(indexij),
					    ilproc=eeList.i2proc(indexij),
					    klproc=Pc.i2proc(indexkl);

					double ijval, ilval, klval;
					
					if(procno==ijproc){
						ijval=eeList.data[indexij-eeList.startpos];
						if(procno!=dataproc){
							MPI_Send(&ijval,1,MPI_DOUBLE,dataproc,0,eeInteract.comm);
						}
						
					}
					else if(procno==dataproc){
						MPI_Recv(&ijval,1,MPI_DOUBLE,ijproc,0,eeInteract.comm,
								MPI_STATUS_IGNORE);
					}
					
					if(procno==ilproc) {
						ilval=eeList.data[indexil-eeList.startpos];
						if(procno!=dataproc){
							MPI_Send(&ilval,1,MPI_DOUBLE,dataproc,0,eeInteract.comm);
						}
					}
					else if(procno==dataproc){
						MPI_Recv(&ilval,1,MPI_DOUBLE,ilproc,0,eeInteract.comm,
								MPI_STATUS_IGNORE);
					}
					if(procno==klproc){
						klval=Pc.data[indexkl-Pc.startpos];
						if(procno!=dataproc){
							MPI_Send(&klval,1,MPI_DOUBLE,dataproc,0,eeInteract.comm);
						}
					}
					else if(procno==dataproc){
						MPI_Recv(&klval,1,MPI_DOUBLE,klproc,0,eeInteract.comm,
								MPI_STATUS_IGNORE);
					}
					
					if(procno==dataproc)sub+=klval*(ijval-0.5*ilval);
					
				}
			}
			if(eeInteract.procno==dataproc) eeInteract.data[eeInteract.rc2i(i,j)-eeInteract.startpos]=sub;
		}
	}
}

DistributedMatrix HFPar::computeNewDensity(DistributedMatrix & C){
	DistributedMatrix newDen(C.rows,C.cols,C.comm);

	int procno= newDen.procno;

	for(int i = 0; i < C.rows; i++){
		for(int j = 0; j < C.cols; j++){
			for(int k = 0; k < (int)(electroncount/2);k++){
				int ijindex=newDen.rc2i(i,j);
				int ijproc=newDen.i2proc(ijindex);
				int ikindex=newDen.rc2i(i,k);
				int ikproc=newDen.i2proc(ikindex);
				int jkindex=newDen.rc2i(j,k);
				int jkproc=newDen.i2proc(jkindex);

				double ikval, jkval;

				if(ikproc==procno){
					ikval=C.data[ikindex-C.startpos];
					if(procno!=ijproc){
						MPI_Send(&ikval,1,MPI_DOUBLE,ijproc,0,C.comm);
					}
				}
				else if(procno==ijproc){
					MPI_Recv(&ikval,1,MPI_DOUBLE,ikproc,0,C.comm,MPI_STATUS_IGNORE);
				}
				if(jkproc==procno){
					jkval=C.data[jkindex-C.startpos];
					if(procno!=ijproc){
						MPI_Send(&jkval,1,MPI_DOUBLE,ijproc,0,C.comm);
					}
				}
				else if(procno==ijproc){
					MPI_Recv(&jkval,1,MPI_DOUBLE,jkproc,0,C.comm,MPI_STATUS_IGNORE);
				}

				if(procno==ijproc){
					newDen.data[ijindex-newDen.startpos]+=2*ikval*jkval;
				}

			}
		}
	}
	return newDen;
}

void HFPar::computeDensityDifference(DistributedMatrix & P_Old, DistributedMatrix & P_New){
	double delta=0.0;

	for(int i = 0; i < P_Old.data.size();i++){
		double diff=P_Old.data[i]-P_New.data[i];
		delta+=diff*diff;
	}
	MPI_Allreduce(&delta,&eDiff,1,MPI_DOUBLE,MPI_SUM,P_Old.comm);

	eDiff=std::sqrt(eDiff);
}

double HFPar::calculateEnergyTotal(Molecule & mol){
	return calculateENuc(mol)+calculateEElectron();
}

double HFPar::calculateENuc(Molecule & mol){
	double En = 0.0;
	for(int i = 0; i < mol.atoms.size(); i++){
		for(int j = i+1; j < mol.atoms.size(); j++){
			Atom & a1=mol.atoms[i], & a2=mol.atoms[j];
			double dprod=0.0;
			for(int i = 0; i < 3; i++){
				double diff=a1.pos[i]-a2.pos[i];
				dprod+=diff*diff;
			}
			double Z1=a1.atnum, Z2=a2.atnum;
			En+=Z1*Z2/(std::sqrt(dprod));
		}
	}
	return En;
}

double HFPar::calculateEElectron(){
	double E = 0.0;
	double lE=0.0;
	for(int i = 0; i < P.data.size();i++){
		lE+=0.5*P.data[i]*(cH.data[i]+F.data[i]);
	}
	MPI_Allreduce(&lE,&E,1,MPI_DOUBLE,MPI_SUM,P.comm);
	return E;
}
