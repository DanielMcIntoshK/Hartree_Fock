#ifndef HF__H_
#define HF__H_
#include "Matrix.h"
#include "Atom.h"
#include "Integrals.h"

class HF{
public:
	HF(Matrix & _S, Matrix & _X, Matrix & _cH, list4D ld);

	double calculateEnergy(Molecule & mol, Matrix & P_init, int charge, bool verbose);

	void HFStep();

	void computeInteraction(Matrix & P_i);

	Matrix computeNewDensity(Matrix & C);
	void computeDensityDifference(Matrix & P_Old, Matrix & P_New);

	double calculateEnergyTotal(Molecule& mol);
private:
	Matrix S, X, Xt, cH,F;

	list4D eeList;

	Matrix P,eeI;

	int dim;

	int electroncount;

	double eDiff;

	double calculateENuc(Molecule & mol);
	double calculateEElectron();

	bool verbose;
};

#endif

