#ifndef ATOM__H_
#define ATOM__H_
#include <vector>
#include <string>

struct pg{
	pg():a{0.0},c{0.0}{}
	pg(double _a, double _c):a{_a},c{_c}{}

	double a;
	double c;

	void generateNorm(std::vector<int> &qnums);

	double norm;
};

struct Orbital{
	void addPG(pg newpg) {newpg.generateNorm(qnums);pgs.push_back(newpg);}
	
	std::vector<pg> pgs; 
	std::vector<int> qnums=std::vector<int>{0,0,0};
	std::vector<double> pos=std::vector<double>{0,0,0};
	unsigned char type;

	Orbital getrel(int off, int dir){
		Orbital no;
		no.pgs=pgs;
		no.pos=pos;
		no.qnums=qnums;
		no.qnums[dir]+=off;
		no.type=type;
		return no;
	}

	void printorbital();
};

class Atom{
public:
	Atom(int _atnum, std::vector<double> ps=std::vector<double>{0,0,0}):atnum{_atnum},pos{ps}{}
	~Atom(){}

	void buildorbitals(std::string basisfile);
	void printorbitals();

	int atnum;
	std::vector<Orbital> orbitals;

	std::vector<double> pos;
	
	std::string getName() {return namelist[atnum-1];}

	static std::vector<std::string> namelist;
private:
	std::vector<Orbital> constructOrbitals(char type);
};

struct Molecule{
	std::vector<Atom> atoms;
	int charge, spinmult;

	int orbitalcount;

	void init(std::string bset);
	Orbital& getOrbital(int n);
};

#endif
