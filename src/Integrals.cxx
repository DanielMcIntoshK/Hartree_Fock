#include "Integrals.h"
#include "Constants.h"
#include "MyMath.h"
#include <vector>
#include <cmath>
#include <iostream>

double overlap(Orbital o1, Orbital o2){
	double Sval=0.0;
	for(int i=0; i< o1.pgs.size();i++){
		for(int j=0; j< o2.pgs.size();j++){
			pg &pg1=o1.pgs[i],&pg2=o2.pgs[j];
			std::vector<double> gprod=gaussianProduct(pg1.a,o1.pos,pg2.a,o2.pos);
			
			//double tval=pg1.c*pg2.c;
			double tval=pg1.c*pg2.c;
			for(int d = 0; d < 3; d++){ 
				double dirov=directionaloverlap(o1,o2,gprod[d],i,j,d);
				tval*=dirov;
			}
		        tval*=gprod[3]*pg1.norm*pg2.norm;
			double normalization=std::pow(MyPI/(pg1.a+pg2.a),3.0/2.0);
			tval*=normalization;
			Sval+=tval;	
		}
	}
	return Sval;
	
}

double kinetic(Orbital o1, Orbital o2){
	double Kval=0;
	for(int i = 0; i < o1.pgs.size();i++){
		for(int j = 0; j < o2.pgs.size();j++){
			pg &pg1=o1.pgs[i], &pg2=o2.pgs[j];

			double Onorm = std::pow(MyPI/(pg1.a+pg2.a),3.0/2.0);

			std::vector<double> gprod=gaussianProduct(pg1.a,o1.pos,pg2.a,o2.pos);
		
			double kij=0;
			for(int d = 0; d < 3; d++){
				kij+=directionalKinetic(o1, o2, gprod[d], i, j, d, gprod,Onorm); 
			}
			kij*=pg1.c*pg2.c*pg1.norm*pg2.norm;

			Kval+=kij;
		}
	}
	return Kval;
}

double nuclear(Orbital o1, Orbital o2, Atom & a){
	double Vn=0;
	double atomcharge = a.atnum;
	for(int pgi = 0; pgi < o1.pgs.size(); pgi++){
	for( int pgj = 0; pgj < o2.pgs.size();pgj++){
		pg &pg1 = o1.pgs[pgi], &pg2 = o2.pgs[pgj];
		double nuc=0.0;

		double g = pg1.a + pg2.a;
		double eps=1.0/(4*g);

		std::vector<double> gprod=gaussianProduct(pg1.a,o1.pos,pg2.a,o2.pos);
		

		for(int l = 0; l < (o1.qnums[0]+o2.qnums[0]+1);l++){
		for(int r = 0; r < (int)(l/2)+1;r++){
		for(int i = 0; i < (int)((l-2*r)/2)+1;i++){
			double Ax=aexpansion(l,r,i,o1.qnums[0],o2.qnums[0],o1.pos[0],o2.pos[0],a.pos[0],gprod[0],eps);
			for(int m = 0; m < (o1.qnums[1]+o2.qnums[1]+1);m++){
			for(int s = 0; s < (int)(m/2)+1;s++){
			for(int j = 0; j < (int)((m-2*s)/2)+1;j++){
				double Ay=aexpansion(m,s,j,o1.qnums[1],o2.qnums[1],o1.pos[1],o2.pos[1],a.pos[1],gprod[1],eps);
				for(int n = 0; n < (o1.qnums[2]+o2.qnums[2]+1);n++){
				for(int t = 0; t < (int)(n/2)+1;t++){
				for(int k = 0; k < (int)((n-2*t)/2)+1; k++){
					double Az=aexpansion(n,t,k,o1.qnums[2],o2.qnums[2],o1.pos[2],o2.pos[2],a.pos[2],gprod[2],eps);

					double nu = l+m+n-2*(r+s+t)-(i+j+k);

					double dprod = 0;
					
					for(int d = 0; d < 3; d++){
						dprod += (gprod[d]-a.pos[d])*(gprod[d]-a.pos[d]);
					}
					double ff = boysFunction(nu,g*dprod);
					nuc+= Ax*Ay*Az*ff;
				}}}
			}}}
		}}}

		nuc*=-atomcharge*pg1.norm*pg2.norm*gprod[3]*2*MyPI/g;
		Vn+=nuc*pg1.c*pg2.c;
	}}
	return Vn;
}

double electronic(Orbital o1, Orbital o2, Orbital o3, Orbital o4){
	double G=0.0;
	for(int pgi = 0; pgi < o1.pgs.size(); pgi++){
	for(int pgj = 0; pgj < o2.pgs.size(); pgj++){
	for(int pgk = 0; pgk < o3.pgs.size(); pgk++){	
	for(int pgl = 0; pgl < o4.pgs.size(); pgl++){
		double sub=0.0;

		pg &pg1=o1.pgs[pgi], &pg2=o2.pgs[pgj], &pg3=o3.pgs[pgk], &pg4=o4.pgs[pgl];
		
		std::vector<int> &a=o1.qnums, &b=o2.qnums, &c=o3.qnums, &d=o4.qnums;

		double g1 = pg1.a + pg2.a, g2 = pg3.a + pg4.a;

		std::vector<double> gprod1 = gaussianProduct(pg1.a,o1.pos,pg2.a,o2.pos);
		std::vector<double> gprod2 = gaussianProduct(pg3.a,o3.pos,pg4.a,o4.pos);

		double delta = 1.0/(4.0*g1) + 1.0/(4.0*g2);
		
		for(int l = 0; l < a[0]+b[0]+1;l++){
		for(int r = 0; r < (int)(l/2)+1;r++){
		for(int ll = 0; ll < c[0]+d[0]+1;ll++){
		for(int rr = 0; rr < (int)(ll/2)+1;rr++){
		for(int i = 0; i < (int)((l+ll-2*r-2*rr)/2)+1;i++){
			double Bx = bexpansion(l,ll,r,rr,i,a[0],b[0],c[0],d[0],o1.pos[0],o2.pos[0],gprod1[0],
					o3.pos[0],o4.pos[0],gprod2[0],g1,g2,delta);
			for(int m = 0; m < a[1]+b[1]+1;m++){
			for(int s = 0; s < (int)(m/2)+1;s++){
			for(int mm = 0; mm < c[1]+d[1]+1;mm++){
			for(int ss = 0; ss < (int)(mm/2)+1;ss++){
			for(int j = 0; j < (int)((m+mm-2*s-2*ss)/2)+1;j++){
				double By = bexpansion(m,mm,s,ss,j,a[1],b[1],c[1],d[1],o1.pos[1],o2.pos[1],gprod1[1],
						o3.pos[1],o4.pos[1],gprod2[1],g1,g2,delta);
				for(int n = 0; n < a[2]+b[2]+1;n++){
				for(int t = 0; t < (int)(n/2)+1;t++){
				for(int nn = 0; nn < c[2]+d[2]+1;nn++){
				for(int tt = 0; tt < (int)(nn/2)+1;tt++){
				for(int k = 0; k < (int)((n+nn-2*t-2*tt)/2)+1;k++){
					double Bz = bexpansion(n,nn,t,tt,k,a[2],b[2],c[2],d[2],o1.pos[2],o2.pos[2],gprod1[2],
							o3.pos[2],o4.pos[2],gprod2[2],g1,g2,delta);

					double nu = l+ll+m+mm+n+nn-2*(r+rr+s+ss+t+tt)-(i+j+k);

					double dprod=0.0;
					for(int i = 0; i < 3; i++) dprod+=(gprod1[i]-gprod2[i])*(gprod1[i]-gprod2[i]);
					dprod/=4.0*delta;

					double ff = boysFunction(nu,dprod);
					
					sub+= Bx*By*Bz*ff;
				}}}}}
				
			}}}}}
		}}}}}
		sub*=pg1.c*pg2.c*pg3.c*pg4.c;
		G+=sub*pg1.norm*pg2.norm*pg3.norm*pg4.norm*gprod1[3]*gprod2[3]*2*std::pow(MyPI,2)/(g1*g2)*std::sqrt(MyPI/(g1+g2));
	}}}}
	return G;
}

double directionaloverlap(Orbital o1, Orbital o2, double center, int pg1, int pg2, int dir){
	double overlap=0.0;
	for(int i = 0; i <= o1.qnums[dir]; i++){
		for(int j = 0; j <= o2.qnums[dir]; j++){
			if(((i+j)%2)!=0)continue;
			double comb1=comb(o1.qnums[dir],i);
			double comb2=comb(o2.qnums[dir],j);
			double facpart=fac2(i+j-1);
			double quotent=std::pow(2*(o1.pgs[pg1].a+o2.pgs[pg2].a),(double)(i+j)/2.0);
			double poly1 = std::pow(center-o1.pos[dir],o1.qnums[dir]-i);
			double poly2 = std::pow(center-o2.pos[dir],o2.qnums[dir]-j);

			overlap+=comb1*comb2*facpart*poly1*poly2/quotent;
		}
	}
	return overlap;
}

double directionalKinetic(Orbital o1, Orbital o2, double center, int pg1, int pg2, int dir, std::vector<double> gprod, double Onorm){
	typedef std::vector<int> iv;
	double recursion=0;
	recursion+=o1.qnums[dir]*o2.qnums[dir]*directionaloverlap(o1.getrel(-1,dir),o2.getrel(-1,dir),center,pg1,pg2,dir);
	recursion+=-2*o1.pgs[pg1].a*o2.qnums[dir]*directionaloverlap(o1.getrel(1,dir),o2.getrel(-1,dir),center,pg1,pg2,dir);
	recursion+=-2*o1.qnums[dir]*o2.pgs[pg2].a*directionaloverlap(o1.getrel(-1,dir),o2.getrel(1,dir),center,pg1,pg2,dir);
	recursion+=4*o1.pgs[pg1].a*o2.pgs[pg2].a*directionaloverlap(o1.getrel(1,dir),o2.getrel(1,dir),center,pg1,pg2,dir);
	recursion*=.5;

	double kv=Onorm*gprod[3]*recursion;

	for(int i = 0; i<3;i++){
		if(i!=dir){
			kv*=directionaloverlap(o1,o2,gprod[i],pg1,pg2,i);
		}
	}

	

	return kv;

}

std::vector<double> gaussianProduct(double a1, std::vector<double> r1, double a2, std::vector<double> r2){
	std::vector<double> r{0.0,0.0,0.0,0.0};
	for(int i = 0; i < 3; i++){r[i]=(a1*r1[i]+a2*r2[i])/(a1+a2);}
	for(int i = 0; i < 3; i++){r[3]+=(r1[i]-r2[i])*(r1[i]-r2[i]);}
	r[3]*= -a1*a2/(a1+a2);
	r[3]=std::exp(r[3]);
	return r;
}

Matrix computeSMatrix(Molecule &mol){
	int dim=mol.orbitalcount;
	Matrix S(dim,dim);
	
	int atcnt=0;
	int ocnt=0;
	for(int r = 0; r < dim; r++){
		for(int c=0; c < dim; c++){
			S.setElement(r,c,overlap(mol.getOrbital(r),mol.getOrbital(c)));
		}	
	}
	return S;
}

Matrix computeKineticMatrix(Molecule & mol){
	int dim = mol.orbitalcount;
	Matrix cH(dim,dim);

	int atcnt=0;
	int ocnt=0;
	for(int r = 0; r < dim; r++){
		for(int c= 0; c < dim; c++){
			cH.setElement(r,c,kinetic(mol.getOrbital(r),mol.getOrbital(c)));
		}
	}
	return cH;
}

Matrix computeNuclearMatrix(Molecule & mol){
	int dim = mol.orbitalcount;
	Matrix N(dim,dim);

	for(int a = 0; a < mol.atoms.size(); a++){
		Matrix AtomMat(dim,dim);
		int atcnt=0;
		int ocnt=0;
		for(int r = 0; r < dim; r++){
			for(int c = 0; c < dim; c++){
				AtomMat.setElement(r,c,nuclear(mol.getOrbital(r),mol.getOrbital(c),mol.atoms[a]));
			}
		}
		std::cout << "NUCLEAR MATRIX " << mol.atoms[a].getName() << ":\n";
		AtomMat.printMatrix();
		N=Matrix::matAdd(N,AtomMat);
	}


	return N;
}

Matrix computeCoreHamiltonianMatrix(Molecule & mol){
	int dim=mol.orbitalcount;
	Matrix Hc(dim,dim), K(dim,dim), N(dim,dim);
	std::cout << "COMPUTEING CORE HAMILTONIAN\n";
	
	K=computeKineticMatrix(mol);

	std::cout << "Kinetic Matrix:\n";
	K.printMatrix();

	N=computeNuclearMatrix(mol);
	
	std::cout << "Nuclear Matrix:\n";
	N.printMatrix();

	Hc=Matrix::matAdd(K,N);

	return Hc;
}

list4D computeEEMatricies(Molecule &mol){
	list4D ld;
	int dim=mol.orbitalcount;
	ld.resize(dim);
	for(int i = 0; i < dim; i++){
		ld[i].resize(dim);
		for(int j = 0; j < dim; j++){
			ld[i][j].resize(dim);
			for(int k = 0; k < dim; k++){
				ld[i][j][k].resize(dim);
			}
		}
	}

	for(int i = 0; i < dim; i++){
	for(int j = 0; j < dim; j++){
	for(int k = 0; k < dim; k++){
	for(int l = 0; l < dim; l++){
		double ev=electronic(mol.getOrbital(i),mol.getOrbital(j), mol.getOrbital(k), mol.getOrbital(l));

		//std::cout << "("<<i+1 <<","<<j+1<<","<<k+1<<","<<l+1<<")"<<"  " << ev << std::endl; 

		ld[i][j][k][l]=ev;
		
	}}}}
	return ld;
}

