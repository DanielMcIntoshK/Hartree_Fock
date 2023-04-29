#include "HF.h"
#include <iostream>
#include <cmath>

HFPar::HFPar(DistributedMatrix & _S, DistributedMatrix & _X, DistributedMatrix & _cH,
		DistributedMatrix & ld,std::vector<double> _eeL):S{_S},X{_X},cH{_cH},eeList{ld},eeListS{_eeL}{
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

	if(S.procno==0&&verbose)std::cout << "ELECTRONS: " << electroncount<<std::endl;

	dim = mol.orbitalcount;

	P=P_init;

	int stepCount = 1;

	do{
		if(S.procno==0) std::cout << std::endl << "SCF STEP " << stepCount++ << std::endl;
		HFStep();
		double totalE = calculateEnergyTotal(mol);
		if(S.procno==0) std::cout << "ENERGY DIFF: " << eDiff << std::endl << "TOTAL ENERGY: " << totalE << std::endl;
	}while(eDiff > 1.0E-12 && stepCount < 10000);
	if(S.procno==0){
	if(stepCount>=10000) std::cout << "\nSCF DID NOT CONVERGE\n";
	else std::cout << "\nSCF CONVERGED\n";
	}
	return calculateEnergyTotal(mol);
}

void HFPar::HFStep(){
	
	//computeInteraction(P);
	computeInteractionStatic(P);
	
	F=DistributedMatrix::matAdd(cH,eeInteract);

	DistributedMatrix Fx = DistributedMatrix::matMul(Xt,F);
	Fx=DistributedMatrix::matMul(Fx,X);

	DistributedEigenSolver es(Fx,0.000001);
	DistributedEigenSolver::EigenData ed = es.calculateEigens();
	
	
	std::vector<int> sl=sortlist(ed.eigenVals);
	ed.eigenVecs.sort(sl);
	
	
	DistributedMatrix C = DistributedMatrix::matMul(X,ed.eigenVecs);
	
	DistributedMatrix P_New=computeNewDensity(C);

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
					int indexil=i*dim*dim*dim+l*dim*dim+k*dim+j;
					
					int indexkl=Pc.rc2i(k,l);

					int ijproc=eeList.i2proc(indexij),
					    ilproc=eeList.i2proc(indexil),
					    klproc=Pc.i2proc(indexkl);

					double ijval=0.0, ilval=0.0, klval=0.0;
					
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
					
					if(procno==dataproc){
						sub+=klval*(ijval-0.5*ilval);
					}
					
				}
			}
			if(eeInteract.procno==dataproc) {
				eeInteract.data[eeInteract.rc2i(i,j)-eeInteract.startpos]=sub;

			}
		}
	}
}

void HFPar::computeInteractionStatic(DistributedMatrix &Pc){
	eeInteract=DistributedMatrix(dim,dim,Pc.comm);
	std::vector<double> pcData=Pc.gatherMat();
	int procno=Pc.procno;
	for(int i = 0; i < dim; i++){
		for(int j = 0; j < dim; j++){
			double sub=0.0;
			for(int k = 0; k < dim; k++){
			for(int l = 0; l < dim; l++){
				long ij=i*dim*dim*dim+j*dim*dim+k*dim+l;
				long il=i*dim*dim*dim+l*dim*dim+k*dim+j;
				int kl=Pc.rc2i(k,l);
				
				double ijval=eeListS[ij], ilval=eeListS[il];
			        double	klval=pcData[kl];

				sub+=klval*(ijval-0.5*ilval);	
			}
			}
			if(procno==Pc.i2proc(Pc.rc2i(i,j))){
				eeInteract.data[eeInteract.rc2i(i,j)-eeInteract.startpos]=sub;
			}
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

	eDiff=std::sqrt(eDiff/4.0);
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
