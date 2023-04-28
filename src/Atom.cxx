#include "Atom.h"
#include "Constants.h"
#include "MyMath.h"
#include <fstream>
#include <iostream>
#include <cmath>

void pg::generateNorm(std::vector<int> &qnums){
	norm=std::pow((2.0*a/MyPI),3.0/4.0);
	norm*=std::pow(4.0*a, (double)(qnums[0]+qnums[1]+qnums[2])/2.0);
	int facstuff=(fac2(2*qnums[0]-1)*fac2(2*qnums[1]-1)*fac2(2*qnums[2]-1));
	norm/=std::sqrt(facstuff);
}

void Orbital::printorbital(){
	std::string names="SPD";
	std::cout << qnums[0] << " " << qnums[1] << " " << qnums[2] << " : " << names[qnums[0]+qnums[1]+qnums[2]]<<" ("<<pos[0]<<","<<pos[1]<<","<<pos[3]<<")\n";
	for(auto &a: pgs){
		std::cout << a.a << " " << a.c << std::endl;
	}
}

void Atom::buildorbitals(std::string basisfile){
	std::ifstream inFile(basisfile);
	if(!inFile){
		std::cout << basisfile << " IS NOT A FILE\n";
	}

	std::string atomname=namelist[atnum-1];	

	std::string linestr;
	bool loading=false;
	int pgcount=-1;
	char orbitaltype;
	std::vector<Orbital> tempOrbitals;

	bool hybrid=false;
	while(!inFile.eof()){
		std::getline(inFile,linestr);
		if(linestr==atomname) {loading=true;continue;}
		if(loading){
			if(linestr.empty()) break;
			if(pgcount==-1){
				orbitaltype=linestr[0];
				hybrid=orbitaltype=='L';
				pgcount=std::stoi(linestr.substr(1))-1;
				tempOrbitals=constructOrbitals(orbitaltype);
			}	
			else{
				
				linestr=linestr.substr(1);
				int pos=linestr.find_first_not_of(' ');
				linestr=linestr.substr(pos);
				pos=linestr.find_first_of(' ');
				if(!hybrid){
					double a=std::stod(linestr.substr(0,pos)),
				    		c=std::stod(linestr.substr(pos));
					for(auto &t:tempOrbitals) t.addPG(pg(a,c));
				}
				else{
					double a=std::stod(linestr.substr(0,pos));
					linestr=linestr.substr(linestr.find_first_not_of(' '));
					linestr=linestr.substr(linestr.find_first_of(' '));
					linestr=linestr.substr(linestr.find_first_not_of(' '));
					pos=linestr.find_first_of(' ');
					double cs=std::stod(linestr.substr(0,pos)),
						cp=std::stod(linestr.substr(pos));
					for(auto &t:tempOrbitals) {
						double ctouse=(t.type==0)?cs:cp;
						t.addPG(pg(a,ctouse));
					}	
				}
				if(--pgcount == -1){
					for(auto a:tempOrbitals) orbitals.push_back(a);	
				}
			}
		}
	}

}

void Atom::printorbitals(){
	std::cout << "ATOM: " << namelist[atnum-1] << std::endl;
	for(auto &a:orbitals){
		a.printorbital();
	}
	std::cout << std::endl;
}

std::vector<Orbital> Atom::constructOrbitals(char type){
	std::vector<Orbital> orbitals;
	switch(type){
	case 'S':{
		Orbital o;
		o.type=0;
		orbitals.push_back(o);
	}break;
	case 'P':{
		for(int i = 0; i < 3; i++){
			Orbital o;
			o.qnums[i]=1;
			o.type=1;
			orbitals.push_back(o);
		}
	}break;
	case 'L':{
		std::vector<Orbital> sorbs=constructOrbitals('S');
		std::vector<Orbital> porbs=constructOrbitals('P');
		orbitals.push_back(sorbs[0]);
		for(auto &o: porbs) orbitals.push_back(o);
	}break;
	case 'D':{
		for(int i = 0; i < 3; i++){
			Orbital o;
			o.qnums[i]=2;
			o.type=2;
			orbitals.push_back(o);
		}
		for(int i = 0; i < 3; i++){
			Orbital o;
			o.qnums=std::vector<int>{1,1,1};
			o.qnums[i]=0;
			o.type=2;
			orbitals.push_back(o);
		}
	}break;
	}
	for(auto &o: orbitals){o.pos=pos;}
	return orbitals;
}

void Molecule::loadMolecule(std::string fileName){
	std::ifstream inFile(fileName);
	std::string input;
	std::getline(inFile, input);
	int atomCount=std::stoi(input);
	std::cout << atomCount << std::endl;
	std::getline(inFile, input);
	for(int i = 0; i < atomCount; i++){
		std::getline(inFile,input);
		std::string at=strip(input),
			px=strip(input),
			py=strip(input),
			pz=strip(input);
		//std::cout << "ATOMINFO:\n" <<"\t"<<at<<"\n\t"<<px <<"\n\t"<<py<<"\n\t"<<pz<<std::endl;
		int atnum;
		for(int i = 0; i <Atom::atabrev.size();i++){
			if(Atom::atabrev[i]==at) {atnum=i+1;break;}
		}	
		std::vector<double> atpos{std::stod(px),std::stod(py),std::stod(pz)};
		atoms.push_back(Atom(atnum,atpos));
	}
}

void Molecule::init(std::string bset){
	orbitalcount=0;
	for(auto & A: atoms){
		A.buildorbitals(bset);
		orbitalcount+=A.orbitals.size();
	}
}

Orbital & Molecule::getOrbital(int n){
	int ocnt=0;
	for(auto & a: atoms){
		for(auto &o:a.orbitals){
			if(ocnt++==n) return o;
		}
	}
	std::cout << "ORBITAL OVERFLOW\n";
	return atoms[0].orbitals[0];
}

void Molecule::printMol(){
	for(auto & a: atoms){
		std::cout << Atom::atabrev[a.atnum-1] << " " << a.pos[0] << " " <<
			a.pos[1] << " " << a.pos[2] << std::endl;
	}
}

std::string Molecule::strip(std::string & str){
	//std::cout << "1:"<< str << std::endl;
	int psf = str.find_first_not_of(' ');
	str=str.substr(psf);
	//std::cout << "2:"<< str << std::endl;
	int pfe = str.find_first_of(' ');
	if(pfe==std::string::npos) return str;
	std::string retstr=str.substr(0,pfe);
	str=str.substr(pfe);
	//std::cout << "3:"<< str << std::endl;
	int end = str.find_first_not_of(' ');
	str=str.substr(end);
	//std::cout << "4:"<<str << std::endl;

	//std::cout <<"R:"<<retstr<<std::endl;
	return retstr;
}

std::vector<std::string> Atom::namelist = std::vector<std::string>{"HYDROGEN","HELIUM","LITHIUM","BERYLLIUM","BORON","CARBON","NITROGEN","OXYGEN","FLUORINE","NEON","SODIUM","MAGNESIUM","ALUMINIUM","SILICON","PHOSPHORUS","SULFUR","CHLORINE","ARGON"};

std::vector<std::string> Atom::atabrev = std::vector<std::string>{"H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar"};
