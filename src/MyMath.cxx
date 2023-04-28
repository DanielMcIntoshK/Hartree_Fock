#include "MyMath.h"
#include <cmath>
#include <iostream>

int fac(int n){
	if(n<=1)return 1;
	return n*fac(n-1);
}

int fac2(int n){
	if(n<=1)return 1;
	return n*fac2(n-2);

}

int comb(int i, int j){
	return fac(i)/(fac(j)*fac(i-j));
}

double binom(double x, double y){
	return tgamma(x+1)/(tgamma(y+1)*tgamma(x-y+1));
}

/*
double gammaIncomplete(double a, double x, double t, bool verbose){
	if(a/x < 0.1){ 
		return tgamma(a);
	}
	double ina=1.0/a;
	double epow=std::pow(x,a);

	double incomplete= ina*epow*hypergeometric(a,1+a,-x,t);

	return incomplete;
}
*/

double gammaIncompleteNew(double p, double x){
	double a, arg, c, f, value;
	double e = 1.0E-09, uflo=1.0E-37;

	if(x<=0.0) {
		std::cout << "GAMMA INCOMPLETE " << p << " " << x << std::endl;
		return 0.0;
	}
	if(p<=0.0){ 
		std::cout << "GAMMA INCOMPLETE " << p << " " << x << std::endl;
		return 0.0;
	}
	
	arg = p * std::log(x)-lgamma(p+1.0)-x;

	if(arg<log(uflo)) return 0.0;

	f=std::exp(arg);

	if(f==0.0) {
		std::cout << "GAMMA INCOMPLETE " << p << " " << x << std::endl;
		
		return 0.0;
	}

	c=1.0;
	value=1.0;
	a=p;

	while(true){
		a=a+1.0;
		c=c*x/a;
		value=value+c;

		if(c<=e*value) break;
	}

	return value*f;
}

double hypergeometric(double a, double b, double z, double t){
	double val = 1;
	double diff = 0.0;

	int k=0;

	double ck=1;

	double ak = 1.0;
	double bk = 1.0;

	double zk=1.0;

	do{
		ak*=a+k;
		bk*=b+k;
		zk*=z;
		
		ck*=++k;
			
		diff=ak*zk/(bk*ck);	

		val+=diff;
	}while(std::abs(diff)>=t&& k<90);

	return val;
}

double boysFunction(double nu, double x){
	double ff;
	if(x < 1E-8) {
		ff = 1.0/(2.0*nu+1.0)-x/(2.0*nu+3.0);
	}
	else {
		ff=0.5/std::pow(x,nu+.5)*tgamma(nu+0.5)*gammaIncompleteNew(nu+0.5,x);
	}
	return ff;

}

double fexpansion(int j, int l, int m, double a, double b){
	double f = 0.0;
	int lowbound = std::max(0,j-m);
	int highbound = std::min(j,l)+1;
	for(int k = lowbound; k < highbound; k++){
		double sub=1;
		sub*=binom(l,k);
		sub*=binom(m,j-k);
		sub*=std::pow(a,l-k);
		sub*=std::pow(b,m+k-j);

		f+=sub;
	}
	return f;
}

double aexpansion(int l, int r, int i, int l1, int l2, double Ra, double Rb, double Rc, double Rp,double eps){
	double A = 1.0;
	if(l%2==1) A*=-1;
	A*=fexpansion(l,l1,l2,Rp-Ra,Rp-Rb);
	if(i%2==1) A*=-1;
	A*=fac(l);
	A*=std::pow(Rp-Rc,l-2*r-2*i);
	A*=std::pow(eps,r+i);
	A/=fac(r);
	A/=fac(i);
	A/=fac(l-2*r-2*i);

	return A;	
}

double thetaexpansion(int l, int l1, int l2, double a, double b, double r, double g){
	double theta=1.0;
	theta*=fexpansion(l,l1,l2,a,b);
	theta*=fac(l);
	theta*=std::pow(g,r-l);
	theta/=fac(r)*fac(l-2*r);

	return theta;
}

double bexpansion(int l, int ll, int r, int rr, int i, int l1, int l2, int l3, int l4,
		double Ra, double Rb, double Rp, double Rc, double Rd, double Rq, double g1, double g2,double delta){
	double b=1.0;
	b*=(l%2==0)?thetaexpansion(l,l1,l2,Rp-Ra,Rp-Rb,r,g1):-thetaexpansion(l,l1,l2,Rp-Ra,Rp-Rb,r,g1);
	b*=thetaexpansion(ll,l3,l4,Rq-Rc,Rq-Rd,rr,g2);
	b*=(i%2==0)?std::pow(2*delta,2*(r+rr)):-std::pow(2*delta,2*(r+rr));
	b*=fac(l+ll-2*r-2*rr);
	b*=std::pow(delta,i)*std::pow(Rp-Rq,l+ll-2*(r+rr+i));

	double denom=1.0;
	denom*=std::pow(4*delta,l+ll)*fac(i);
	denom*=fac(l+ll-2*(r+rr+i));

	return b/denom;
}




