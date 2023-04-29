#ifndef HF__H_
#define HF__H_
#include "Matrix.h"
#include "Atom.h"
#include "Integrals.h"
#include "DistributedMatrix.h"

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

class HFPar{
public:
	HFPar(DistributedMatrix & _S, DistributedMatrix & _X, DistributedMatrix & _cH, DistributedMatrix & ld,
			std::vector<double> _eeListS);

	double calculateEnergy(Molecule & mol, DistributedMatrix & P_init, int charge, bool verbose);

	void HFStep();

	void computeInteraction(DistributedMatrix & Pc);
	void computeInteractionStatic(DistributedMatrix &Pc);

	DistributedMatrix computeNewDensity(DistributedMatrix & C);

	void computeDensityDifference(DistributedMatrix & P_Old, DistributedMatrix & P_New);

	double calculateEnergyTotal(Molecule & mol);

private:
	DistributedMatrix S, X, Xt, cH, F, P, eeList, eeInteract;

	std::vector<double> eeListS;

	int dim;
	int electroncount;

	double eDiff;
	double calculateENuc(Molecule & mol);
	double calculateEElectron();

	bool verbose;
};

#endif

