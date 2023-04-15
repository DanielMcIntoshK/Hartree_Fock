#include "Atom.h"
#include "Matrix.h"
#include <iostream>

int main(int argc, char ** argv){
	Atom H(1),He(2),Li(3),Be(4),B(5),C(6);
	
	H.buildorbitals("basissets/def2-svp.1.orca");
	He.buildorbitals("basissets/def2-svp.1.orca");
	Li.buildorbitals("basissets/def2-svp.1.orca");
	Be.buildorbitals("basissets/def2-svp.1.orca");
	B.buildorbitals("basissets/def2-svp.1.orca");
	C.buildorbitals("basissets/def2-svp.1.orca");

	H.printorbitals();
	He.printorbitals();
	Li.printorbitals();
	Be.printorbitals();
	B.printorbitals();
	C.printorbitals();

	Molecule H2O;
	H2O.atoms=std::vector<Atom>{Atom(1,{1.4305507,0,0}),Atom(1,{-1.4305507,0,0}),Atom(8,{0,1.1072514})};
	H2O.init("basissets/def2-svp.1.orca");
	std::cout << "H2O: " << H2O.orbitalcount << std::endl;
	Matrix m1(4,3),m2(3,6);

	m1.setData(std::vector<double>{
		1.0, 5.0, 2.9,
		3.4, 2.0, 8.7,
		1.1, 10.3, 14.2,
		0.0, 9.6, 3.0});
	m2.setData(std::vector<double>{
		1.1,9.5,2.2,4.1,7.7,8.8,
		0.3,1.2,1.3,6.6,11.0,12.0,
		1.0,2.0,3.0,4.0,5.0,6.0});

	Matrix m3=Matrix::matMul(m1,m2);

	for(int i = 0; i < m3.rows;i++){
		for(int j=0; j < m3.cols;j++){
			std::cout << m3.data[i*m3.cols+j] << "\t";
		}
		std::cout << std::endl;
	}	

	return 0;
}
