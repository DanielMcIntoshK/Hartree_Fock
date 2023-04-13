#include "Atom.h"
#include <fstream>
#include <iostream>

void Orbital::printorbital(){
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

	while(!inFile.eof()){
		std::getline(inFile,linestr);
		if(linestr==atomname) {loading=true;continue;}
		if(loading){
			if(linestr.empty()) break;
			if(pgcount==-1){
				orbitaltype=linestr[0];
				pgcount=std::stoi(linestr.substr(1))-1;
				orbitals.push_back(Orbital());	
			}	
			else{
				std::cout << linestr<< std::endl;
				linestr=linestr.substr(1);
				int pos=linestr.find_first_not_of(' ');
				linestr=linestr.substr(pos);
				pos=linestr.find_first_of(' ');
				double a=std::stod(linestr.substr(0,pos)),
				    c=std::stod(linestr.substr(pos));
				orbitals.back().addPG(pg(a,c));
				pgcount--;
			}
		}
	}

}

void Atom::printorbitals(){
	std::cout << "ATOM: " << namelist[atnum-1] << std::endl;
	for(auto &a:orbitals){
		a.printorbital();
		std::cout << std::endl;
	}
	std::cout << std::endl;
}


std::vector<std::string> Atom::namelist = std::vector<std::string>{"HYDROGEN","HELIUM","LITHIUM","BERYLLIUM","BORON","CARBON","NITROGEN","OXYGEN","FLUORINE","NEON","SODIUM","MAGNESIUM","ALUMINIUM","SILICON","PHOSPHORUS","SULFUR","CHLORINE","ARGON"};
