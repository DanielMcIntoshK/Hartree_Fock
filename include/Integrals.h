#ifndef INTEGRALS__H_
#define INTEGRALS__H_
#include "Atom.h"
#include "Matrix.h"
#include "DistributedMatrix.h"

typedef std::vector<std::vector<std::vector<std::vector<double>>>> list4D;

double overlap(Orbital o1, Orbital o2);

double kinetic(Orbital o1, Orbital o2);

double nuclear(Orbital o1, Orbital o2, Atom & a);

double electronic(Orbital o1, Orbital o2,Orbital o3, Orbital o4);

double directionaloverlap(Orbital o1, Orbital o2, double center, int pg1, int pg2, int dir);

double directionalKinetic(Orbital o1, Orbital o2, double center, int pg1, int pg2, int dir,std::vector<double> gprod, double Onorm);

std::vector<double> gaussianProduct(double a1, std::vector<double> r1, double a2, std::vector<double> r2);

Matrix computeSMatrix(Molecule &mol);

Matrix computeKineticMatrix(Molecule &mol);

Matrix computeNuclearMatrix(Molecule &mol);

Matrix computeCoreHamiltonianMatrix(Molecule &mol);

list4D computeEEMatricies(Molecule &mol);

DistributedMatrix computeSMatrixPar(Molecule &mol);

DistributedMatrix computeKineticMatrixPar(Molecule & mol);

DistributedMatrix computeNuclearMatrixPar(Molecule & mol);

DistributedMatrix computeCoreHamiltonianMatrixPar(Molecule & mol);

DistributedMatrix computeEEMatriciesPar(Molecule&mol);

#endif

