#include "Atom.h"
#include <fstream>
#include <iostream>

void Orbital::printorbital(){
	std::string names="SPD";
	std::cout << qnums[0] << " " << qnums[1] << " " << qnums[2] << " : " << names[qnums[0]+qnums[1]+qnums[2]]<<std::endl;
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

	while(!inFile.eof()){
		std::getline(inFile,linestr);
		if(linestr==atomname) {loading=true;continue;}
		if(loading){
			if(linestr.empty()) break;
			if(pgcount==-1){
				orbitaltype=linestr[0];
				pgcount=std::stoi(linestr.substr(1))-1;
				tempOrbitals=constructOrbitals(orbitaltype);
			}	
			else{
				std::cout << linestr<< std::endl;
				linestr=linestr.substr(1);
				int pos=linestr.find_first_not_of(' ');
				linestr=linestr.substr(pos);
				pos=linestr.find_first_of(' ');
				double a=std::stod(linestr.substr(0,pos)),
				    c=std::stod(linestr.substr(pos));
				for(auto &t:tempOrbitals) t.addPG(pg(a,c));
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
		orbitals.push_back(o);
	}break;
	case 'P':{
		for(int i = 0; i < 3; i++){
			Orbital o;
			o.qnums[i]=1;
			orbitals.push_back(o);
		}
	}break;
	case 'D':{
		for(int i = 0; i < 3; i++){
			Orbital o;
			o.qnums[i]=2;
			orbitals.push_back(o);
		}
		for(int i = 0; i < 3; i++){
			Orbital o;
			o.qnums=std::vector<int>{1,1,1};
			o.qnums[i]=0;
			orbitals.push_back(o);
		}
	}break;
	}
	return orbitals;
}

void Molecule::init(std::string bset){
	orbitalcount=0;
	for(auto & A: atoms){
		A.buildorbitals(bset);
		orbitalcount+=A.orbitals.size();
	}
}

std::vector<std::string> Atom::namelist = std::vector<std::string>{"HYDROGEN","HELIUM","LITHIUM","BERYLLIUM","BORON","CARBON","NITROGEN","OXYGEN","FLUORINE","NEON","SODIUM","MAGNESIUM","ALUMINIUM","SILICON","PHOSPHORUS","SULFUR","CHLORINE","ARGON"};
