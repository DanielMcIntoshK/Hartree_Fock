#include "DistributedMatrix.h"
#include <set>
#include <iostream>
#include <cmath>

DistributedMatrix::DistributedMatrix(MPI_Comm cm):comm{cm},rows{0},cols{0},size{0}{
	MPI_Comm_size(comm, &nprocs);
	MPI_Comm_rank(comm, &procno);

	blocksizebase=overflow=blocksize=startpos=0;
	//int blocksize = (int)(size/nprocs)+((procno<(size%nprocs))?1:0);
}

DistributedMatrix::DistributedMatrix(int r, int c, MPI_Comm cm):comm{cm},rows{r},cols{c},size{r*c}{
	MPI_Comm_size(comm, &nprocs);
	MPI_Comm_rank(comm, &procno);

	//if(procno>size){procno=-1;return DistributedMatrix(comm);}

	blocksizebase=size/nprocs;
	overflow = size%nprocs;

	blocksize = blocksizebase+((procno<overflow)?1:0);

	startpos=(int)(blocksizebase)*procno+std::min(procno,overflow);

	data.resize(blocksize);
	for(auto &a:data)a=0.0;

}

int DistributedMatrix::i2proc(int i){
	int proc = i/blocksizebase;
	int procextra=i%blocksizebase;
	if(procextra < std::min(overflow,proc)) proc--;
	return proc;
		
}

DistributedMatrix DistributedMatrix::matMul(DistributedMatrix &m1, DistributedMatrix & m2){
	if(m1.procno == -1 || m2.procno==-1){return DistributedMatrix(m1.comm);}
	if(m1.cols!=m2.rows){return DistributedMatrix(m1.comm);}

	DistributedMatrix mm(m1.rows, m2.cols,m1.comm);

	std::vector<MPI_Comm> rowcomms, colcomms;
	rowcomms.resize(m1.rows); colcomms.resize(m2.cols);

	std::set<int> rowinfo,colinfo;

	for(int i = 0; i < m1.data.size();i++){
		rowinfo.insert(m1.i2r(i+m1.startpos));
	}
	for(int i = 0; i < m1.rows; i++){
		int color=(rowinfo.find(i)!=rowinfo.end())?0:1;
		MPI_Comm_split(m1.comm,color,m1.procno,&(rowcomms[i]));
	}
	for(int i = 0; i < m2.data.size();i++){
		colinfo.insert(m2.i2c(i+m2.startpos));
	}
	for(int i = 0; i < m2.cols; i++){
		int color=(colinfo.find(i)!=colinfo.end())?0:1;
		MPI_Comm_split(m2.comm,color,m2.procno,&(colcomms[i]));
	}

	std::vector<double> rmul, cmul, sendvec;
	int sendcount=0;
	std::vector<int> recvcount;
	std::vector<int> displace;
	rmul.resize(m1.cols); cmul.resize(m2.rows);
	int procloadcount=0;
	int dtcount=0;

	for(int r = 0; r < mm.rows; r++){
		for(int c = 0; c < mm.cols; c++){
			int index = mm.rc2i(r,c);
			int procmul = mm.i2proc(index);

			if(rowinfo.find(r)!=rowinfo.end()){
				int recvcur=-1;
				sendcount=0;
				recvcount.clear();
				sendvec.clear();
				displace.clear();
				//gather row
				int gathercount=0;
				for(int i = 0; i < m1.cols; i++){
					//int li=i+m1.startpos;
					int li=m1.rc2i(r,i);
					int lproc=m1.i2proc(li);
					if(lproc!=recvcur){
						if(recvcur==-1)displace.push_back(0);
						else displace.push_back(displace.back()+recvcount.back());
						//displace.push_back(0);
						recvcount.push_back(0);
						recvcur=lproc;
					}
					recvcount.back()++;
					if(lproc==m1.procno){
						int sendpos=li-m1.startpos;
						sendvec.push_back(m1.data[sendpos]);
						sendcount++;
					}
					else{gathercount++;}
				}
				if(gathercount ==0){
					rmul=sendvec;
				}
				else{
					MPI_Allgatherv(sendvec.data(),sendcount,MPI_DOUBLE,rmul.data(),recvcount.data(),
						displace.data(),MPI_DOUBLE,rowcomms[r]);
				}
			}
			if(colinfo.find(c)!=colinfo.end()){
				int recvcur=-1;
				sendcount=0;
				recvcount.clear();
				sendvec.clear();
				displace.clear();
				int gathercount=0;
				for(int i = 0; i < m2.rows; i++){
					int li=m2.rc2i(i,c);
					int lproc=m2.i2proc(li);
					if(lproc!=recvcur){
						if(recvcur==-1)displace.push_back(0);
						else displace.push_back(displace.back()+recvcount.back());
						recvcount.push_back(0);
						recvcur=lproc;
					}
					recvcount.back()++;
					if(lproc==m2.procno){
						int sendpos=li-m2.startpos;
						sendvec.push_back(m2.data[sendpos]);
						sendcount++;
					}
					else gathercount++;
				}
				if(gathercount==0){
					cmul=sendvec;
				}
				else{
					MPI_Allgatherv(sendvec.data(),sendcount,MPI_DOUBLE,cmul.data(),recvcount.data(),
						displace.data(),MPI_DOUBLE,colcomms[c]);
				}
			}
			
			if(procmul==mm.procno){
				double dprod=0.0;
				for(int i = 0; i < m1.cols;i++){
					dprod+=rmul[i]*cmul[i];
				}

				mm.data[dtcount++]=dprod;
			}
			MPI_Barrier(m1.comm);
		}
	}
	return mm;	
}

DistributedMatrix DistributedMatrix::diag(std::vector<double> eigenvals){
	int dim = eigenvals.size();
	DistributedMatrix dm(dim,dim,MPI_COMM_WORLD);
	
	if(dm.procno==-1) return dm;

	for(int i = 0; i < dm.data.size(); i++){
		int index=i+dm.startpos;
		int r=dm.i2r(index);
		int c=dm.i2c(index);
		if(r==c) dm.data[i]=eigenvals[r];
		else dm.data[i]=0.0;
	}
	return dm;
}

DistributedMatrix DistributedMatrix::matAdd(DistributedMatrix &m1, DistributedMatrix & m2){
	if(m1.procno == -1 || m2.procno==-1){return DistributedMatrix(m1.comm);}
	if(m1.rows!=m2.rows || m1.cols!=m2.cols){return DistributedMatrix(m1.comm);}
	DistributedMatrix am(m1.rows,m2.cols,m1.comm);

	for(int i = 0; i < am.data.size();i++){
		am.data[i]=m1.data[i]+m2.data[i];
	}
	return am;
}

DistributedMatrix DistributedMatrix::transpose(){
	DistributedMatrix tm(rows,cols,comm);
	for(int i = 0; i < tm.rows; i++){
		for(int j = i; j < tm.cols; j++){
			int index=rc2i(i,j), tindex=rc2i(j,i);
			int dtproc=i2proc(index), tdtproc=i2proc(tindex);
			if(i==j&&procno==dtproc) tm.data[index-startpos]=data[index-startpos];
			else if(tm.procno==dtproc && tm.procno==tdtproc){
				tm.data[index-startpos]=data[tindex-startpos];
				tm.data[tindex-startpos]=data[index-startpos];
			}
			else if(tm.procno==dtproc){
				double holder = 0.0;
				MPI_Sendrecv(&data[index-startpos], 1, MPI_DOUBLE, tdtproc,0,&holder,
						1,MPI_DOUBLE,tdtproc,0,tm.comm,MPI_STATUS_IGNORE);
				tm.data[index-startpos]=holder;
			}
			else if(tm.procno==tdtproc){
				double holder = 0.0;
				MPI_Sendrecv(&data[tindex-startpos], 1, MPI_DOUBLE, dtproc,0,&holder,
						1,MPI_DOUBLE,dtproc,0,tm.comm,MPI_STATUS_IGNORE);
				tm.data[tindex-startpos]=holder;
			}

		}
	}
	return tm;
}

void DistributedMatrix::identity(){
	for(int i = 0; i < data.size();i++){
		int gpos=i+startpos;
		int gr=i2r(gpos), gc=i2c(gpos);

		if(gr==gc) data[i]=1.0;
		else data[i]=0.0;
	}
}

void DistributedMatrix::printMatrix(){
	std::vector<double> fullMat;
	if(procno==0) fullMat.resize(size);

	std::vector<int> recvcount, displace;
	recvcount.resize(nprocs);
	displace.resize(nprocs);

	for(auto &a: recvcount)a=0;
	for(auto &a: displace)a=0;

	int runningcount=0;
	int currentproc=0;

	for(int i = 0; i < rows; i++){
		for(int j = 0; j < cols; j++){
			int lproc=i2proc(rc2i(i,j));
			recvcount[lproc]++;
			if(currentproc!=lproc){
				displace[lproc]=runningcount;
				currentproc=lproc;
			}
			runningcount++;
		}
	}

	if(procno==0){
	std::cout << "RECVCOUNTS: ";
	for(auto & a: recvcount) std::cout << a << " ";
	std::cout << std::endl;
	
	std::cout << "DISPLACEMENTS: ";
	for(auto & a: displace) std::cout << a << " ";
	std::cout << std::endl;
	}

	MPI_Gatherv(data.data(),blocksize,MPI_DOUBLE,fullMat.data(),recvcount.data(),displace.data(),
			MPI_DOUBLE,0,comm);

	if(procno!=0)return;


	for(int i = 0; i < fullMat.size(); i++){
		std::cout << fullMat[i] << " ";
		if(i%cols==cols-1)std::cout << std::endl;
	}

}


DistributedEigenSolver::DistributedEigenSolver(DistributedMatrix & inMat, double _t):threshold{_t},oMat{inMat},
comm{inMat.comm}
{

}

DistributedEigenSolver::EigenData DistributedEigenSolver::calculateEigens(){
	DistributedMatrix cMat=oMat;
	EigenData ed(oMat.rows,oMat.comm);


	ed.eigenVecs.identity();
	for(int i = 0; i < ed.eigenVals.size(); i++) ed.eigenVals[i]=0.0;

	int citer=0;
	int maxiter=cMat.rows*cMat.rows*cMat.rows;

	findMaxUpperTriangle(cMat);
	do{
		if(cMat.procno==0){
			std::cout << citer << std::endl;
		}
		jacobiRotate(cMat, ed);
		findMaxUpperTriangle(cMat);
		citer++;
	}while(std::fabs(maxvl.val)>threshold && citer < maxiter);
//cMat.printMatrix();
	//std::cout << "GOT OUT " << maxvl.val << " " << threshold  << std::endl;
	if(citer==maxiter&&oMat.procno==0) std::cout << "TIMEOUT\n";
	for(int n = 0; n < cMat.rows; n++){
		int nnproc=cMat.i2proc(cMat.rc2i(n,n));
		
		double diagval=0.0;
		if(cMat.procno==nnproc) diagval=cMat.data[cMat.rc2i(n,n)-cMat.startpos];
		MPI_Bcast(&diagval,1,MPI_DOUBLE,nnproc,cMat.comm);

		ed.eigenVals[n]=diagval;
	}
	ed.eigenVecs=ed.eigenVecs.transpose();
	return ed;
}

void DistributedEigenSolver::findMaxUpperTriangle(DistributedMatrix &cMat){
	valloc lvl;
	lvl.val=0.0;
	lvl.loc=-1;
	for(int i = 0; i < cMat.data.size(); i++){
		int index = i + cMat.startpos;
		int gr=cMat.i2r(index);
		int gc=cMat.i2c(index);
	
		if(gc>gr){
			if(std::fabs(cMat.data[i]) > lvl.val){
				lvl.val=std::fabs(cMat.data[i]);
				lvl.loc = index;
			}	
		}
	}
	MPI_Allreduce(&lvl,&maxvl,1,MPI_DOUBLE_INT,MPI_MAXLOC,cMat.comm);
	
	maxr=cMat.i2r(maxvl.loc);
	maxc=cMat.i2c(maxvl.loc);

	int diagrproc=cMat.i2proc(cMat.rc2i(maxr,maxr)), diagcproc=cMat.i2proc(cMat.rc2i(maxc,maxc)),
		maxtransproc=cMat.i2proc(cMat.rc2i(maxc,maxr));    
	if(diagrproc==cMat.procno){
		diagr=cMat.data[cMat.rc2i(maxr,maxr)-cMat.startpos];
	}
	MPI_Bcast(&diagr,1,MPI_DOUBLE,diagrproc,cMat.comm);
	if(diagcproc==cMat.procno){
		diagc=cMat.data[cMat.rc2i(maxc,maxc)-cMat.startpos];
	}
	MPI_Bcast(&diagc,1,MPI_DOUBLE,diagcproc,cMat.comm);
	if(maxtransproc==cMat.procno){
		maxvalTrans=cMat.data[cMat.rc2i(maxc,maxr)-cMat.startpos];
	}
	MPI_Bcast(&maxvalTrans,1,MPI_DOUBLE,maxtransproc,cMat.comm);
	if(cMat.procno==cMat.i2proc(maxvl.loc)) maxvl.val=cMat.data[maxvl.loc-cMat.startpos];
	MPI_Bcast(&maxvl.val,1,MPI_DOUBLE,cMat.i2proc(maxvl.loc),cMat.comm);

}

void DistributedEigenSolver::jacobiRotate(DistributedMatrix & cMat, EigenData & ed){
	double c, s;
	
	//if(cMat.procno==0)std::cout << maxr << " " << maxc << " " <<maxvl.val << std::endl;
	if(maxvl.val != 0.0){
		double tau, t;
		tau = (diagc-diagr)/(2.0*maxvalTrans);
		if(tau > 0.0000001){
			t = 1.0/(tau+std::sqrt(1.0+tau*tau));
		}
		else{
			t= -1.0/(-tau+std::sqrt(1.0+tau*tau));
		}
		c=1.0/std::sqrt(1+t*t);
		s=c*t;
	}
	else{
		c=1.0;
		s=0.0;
	}
	//std::cout << c << " " << s << diagc << " " << diagr << " " << maxvl.val << std::endl;
	for(int j = 0; j < cMat.rows; j++){
		double iRow, jRow;
		int iind=cMat.rc2i(maxr,j),
		    jind=cMat.rc2i(maxc,j);
		int iproc=cMat.i2proc(iind),
		    jproc=cMat.i2proc(jind);

		//if(cMat.procno==0)std::cout << j << " " <<iproc << " " << jproc<< std::endl;
		if(iproc==cMat.procno) iRow=ed.eigenVecs.data[iind-cMat.startpos];
		if(jproc==cMat.procno) jRow=ed.eigenVecs.data[jind-cMat.startpos];
		

		if(iproc==cMat.procno && jproc != cMat.procno){
			//std::cout << "SENDRECV: " << iproc << " TO " << jproc << std::endl;
			MPI_Sendrecv(&iRow,1,MPI_DOUBLE,jproc,0,&jRow,1,MPI_DOUBLE,jproc,0,
					cMat.comm,MPI_STATUS_IGNORE);
		}
		if(jproc==cMat.procno&& iproc != cMat.procno){
			//std::cout << "SENDRECV: " << jproc << " TO " << iproc << std::endl;
			MPI_Sendrecv(&jRow,1,MPI_DOUBLE,iproc,0,&iRow,1,MPI_DOUBLE,iproc,0,
					cMat.comm,MPI_STATUS_IGNORE);
		}
		if(iproc==cMat.procno){
			ed.eigenVecs.data[iind-cMat.startpos]=iRow*c-jRow*s;
		}
		if(jproc==cMat.procno){
			ed.eigenVecs.data[jind-cMat.startpos]=iRow*s+jRow*c;
		}
	}
	
	//cMat.printMatrix();

	//if(cMat.procno==0)std::cout << "EIGENS DONE\n";
	int diagrproc=cMat.i2proc(cMat.rc2i(maxr,maxr)), diagcproc=cMat.i2proc(cMat.rc2i(maxc,maxc)),
	    maxvlproc=cMat.i2proc(maxvl.loc), maxvltproc=cMat.i2proc(cMat.rc2i(maxc,maxr));

	if(diagrproc==cMat.procno) cMat.data[cMat.rc2i(maxr,maxr)-cMat.startpos]=c*c*diagr-2.0*s*c*maxvl.val+s*s*diagc;
	if(diagcproc==cMat.procno) cMat.data[cMat.rc2i(maxc,maxc)-cMat.startpos]=s*s*diagr+2.0*s*c*maxvl.val+c*c*diagc;
	if(maxvlproc==cMat.procno) cMat.data[cMat.rc2i(maxr,maxc)-cMat.startpos]=0.0;
	if(maxvltproc==cMat.procno) cMat.data[cMat.rc2i(maxc,maxr)-cMat.startpos]=0.0;
	
	
	for(int l = 0; l < cMat.rows; l++){
		if(l!=maxr && l != maxc){
			int rlproc=cMat.i2proc(cMat.rc2i(maxr,l)), lrproc=cMat.i2proc(cMat.rc2i(l,maxr)),
			    clproc=cMat.i2proc(cMat.rc2i(maxc,l)), lcproc=cMat.i2proc(cMat.rc2i(l,maxc));
			
			double eigenil, eigenjk;
			if(cMat.procno==rlproc) eigenil=cMat.data[cMat.rc2i(maxr,l)-cMat.startpos];
			if(cMat.procno==clproc) eigenjk=cMat.data[cMat.rc2i(maxc,l)-cMat.startpos];


			if(cMat.procno==rlproc && cMat.procno !=clproc){
				MPI_Sendrecv(&eigenil,1,MPI_DOUBLE,clproc,0,&eigenjk,1,MPI_DOUBLE,clproc,0,cMat.comm,
						MPI_STATUS_IGNORE);
			}
			if(cMat.procno==clproc && cMat.procno !=rlproc){
				MPI_Sendrecv(&eigenjk,1,MPI_DOUBLE,rlproc,0,&eigenil,1,MPI_DOUBLE,rlproc,0,cMat.comm,
						MPI_STATUS_IGNORE);
			}


			if(cMat.procno==rlproc){
				int index = cMat.rc2i(maxr,l)-cMat.startpos;
				cMat.data[index]=c*eigenil-s*eigenjk;
				if(cMat.procno==lrproc) cMat.data[cMat.rc2i(l,maxr)-cMat.startpos]=cMat.data[index];
				else{
					MPI_Send(&cMat.data[index],1,MPI_DOUBLE,lrproc,0,cMat.comm);
				}
			}
			else if(cMat.procno == lrproc){
				MPI_Recv(&cMat.data[cMat.rc2i(l,maxr)-cMat.startpos],1,MPI_DOUBLE,rlproc,0,cMat.comm,
						MPI_STATUS_IGNORE);
			}

			if(cMat.procno==clproc){
				int index = cMat.rc2i(maxc,l)-cMat.startpos;
				cMat.data[index]=s*eigenil+c*eigenjk;
				if(cMat.procno==lcproc) cMat.data[cMat.rc2i(l,maxc)-cMat.startpos]=cMat.data[index];
				else{
					MPI_Send(&cMat.data[index],1,MPI_DOUBLE,lcproc,0,cMat.comm);
				}
			}
			else if(cMat.procno == lcproc){
				MPI_Recv(&cMat.data[cMat.rc2i(l,maxc)-cMat.startpos],1,MPI_DOUBLE,clproc,0,cMat.comm,
						MPI_STATUS_IGNORE);
			}

		}
	}
	//cMat.printMatrix();
}

