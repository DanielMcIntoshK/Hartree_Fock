#include "Atom.h"

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

	return 0;
}
