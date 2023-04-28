#ifndef MyMath__H_
#define MyMath__H_

int fac(int n);

int fac2(int n);

int comb(int i, int j);

double binom(double x, double y);

//double gammaIncomplete(double a, double x,double t,bool v=false);

double gammaIncompleteNew(double x, double p);

double hypergeometric(double a, double b, double z, double t);

double boysFunction(double nu, double x);

double fexpansion(int j, int l, int m, double a, double b);

double aexpansion(int l, int r, int i, int l1, int l2, double Ra, double Rb, double Rc, double Rp, double eps);

double thetaexpansion(int l, int l1, int l2, double a, double b, double r, double g);

double bexpansion(int l, int ll, int r, int rr, int i, int l1, int l2,int l3, int l4,
		double Ra, double Rb, double Rp, double Rc, double Rd, double Rq, double g1, double g2, double delta);

#endif

