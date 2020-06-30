#include "Modules/Module_Compound.h"
#include "Databases/NuclearMass.h"
#include "CrossSection.h"
#include "CompoundStates.h"

#include <iostream>
#include <iomanip> 
#include <string>
#include <fstream>
#include <sstream>
#include <chrono>
#include <cmath>




namespace Module_Compound{

    void Go(int argc, char *argv[]){
        auto start = std::chrono::steady_clock::now();
        std::cout << std::endl << "Module: compound" << std::endl;
        std::cout << std::endl;
        if(argc < 4){
            std::cout << std::endl << "Not enough input parameters ..." << std::endl;
            PrintHelp();
            exit(1);
        }

        std::string reaction(argv[2]);

        if(argc == 4){
            CalcWidths(reaction, atof(argv[3]));
        }

        if(argc == 5){
            std::string state(argv[4]);
            CalcWidths(reaction, atof(argv[3]), state);
        }

        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << std::endl << "Total calculation time: " << elapsed_seconds.count() << "s\n";

    }

    void CalcWidths(std::string reactionString, double energy, std::string state){
            
        int A = massNumberIntFromString(reactionString);
        int Z = atomicNumberIntFromString(reactionString);
        int pType = pTypeIntFromString(reactionString);
        int parity = parityIntFromString(state);
        double spin = spinDoubleFromString(state);

        std::string energyFile;
        if(!fexists(energyFile.c_str())){
            std::ofstream fp;
            fp.open("output/energyFile.dat");
            fp << energy << std::endl;
            fp.close();
        }

        energyFile="output/energyFile.dat";

        bool forRates = false;
        int entranceState = 0;
        std::vector<int> exitStates(4,-1);

        std::cout <<  "Calculating Compound states for A: " << A << " Z: " << Z << " Projectile: " << pType 
        << " State: " << state << " Spin: " << spin << " Parity: " << parity << std::endl << std::endl;

        CrossSection* xs = new CrossSection(Z,A,pType,energyFile,forRates,entranceState,exitStates);

        if(xs->IsValid())
        {
            CompoundStates compound = xs->CalcCompoundWidth(spin,parity);
            PrintCompoundStates(compound.states);
        }

    }

    void CalcWidths(std::string reactionString, double energy){        
        int A = massNumberIntFromString(reactionString);
        int Z = atomicNumberIntFromString(reactionString);
        int pType = pTypeIntFromString(reactionString);
        std::string energyFile;
        if(!fexists(energyFile.c_str())){
            std::ofstream fp;
            fp.open("output/energyFile.dat");
            fp << energy << std::endl;
            fp.close();
        }

        energyFile="output/energyFile.dat";

        bool forRates = false;
        int entranceState = 0;
        std::vector<int> exitStates(4,-1);

        std::cout <<  "Calculating Compound states for A: " << A << " Z: " << Z << " Projectile: " << pType << std::endl << std::endl;

        CrossSection* xs = new CrossSection(Z,A,pType,energyFile,forRates,entranceState,exitStates);

        if(xs->IsValid())
        {
            CompoundStates compound = xs->CalcCompoundWidth();
            PrintCompoundStates(compound.states);
        }
    }

    void PrintCompoundStates(std::vector<CompoundState> states){
        std::cout     << std::setw(6) <<  std::fixed << "Mass" << "\t"
                      << std::setw(6) <<  std::fixed << "Charge" << "\t"
                      << std::setw(6) <<  std::fixed << "Spin" << "\t"
                      << std::setw(6) <<  std::fixed << "Parity" << "\t"
                      << std::setw(6) <<  std::fixed << "Energy" << "\t"
                      << std::setw(9) <<  std::scientific << "Density" << "\t"
                      << std::setw(6) <<  std::scientific << "g GS Width  " << "\t"
                      << std::setw(6) <<  std::scientific << "n GS Width  " << "\t"
                      << std::setw(6) <<  std::scientific << "p GS Width  " << "\t"
                      << std::setw(6) <<  std::scientific << "a GS Width  " << "\t"
                      << std::setw(6) <<  std::scientific << "g TotalWidth" << "\t"
                      << std::setw(6) <<  std::scientific << "n TotalWidth" << "\t"
                      << std::setw(6) <<  std::scientific << "p TotalWidth" << "\t"
                      << std::setw(6) <<  std::scientific << "a TotalWidth" << "\t"                                            
                      << std::setw(6) <<  std::scientific << "TotalWidth" << "\t"
                      << std::fixed << std::endl;

        for(auto it = states.begin(); it != states.end(); ++it){
            std::cout << std::setw(6) <<  std::fixed << it->A() << "\t"
                      << std::setw(6) <<  std::fixed << it->Z() << "\t"
                      << std::setw(6) <<  std::fixed << it->J() << "\t"
                      << std::setw(6) <<  std::fixed << ((it->P()>0)? "+" :"-" ) << "\t"
                      << std::setw(6) <<  std::fixed << it->E() << "\t"
                      << std::setw(9) <<  std::scientific << it->Density() << "\t"
                      << std::setw(6) <<  std::scientific << it->GroundStateWidth("gamma") << "\t"
                      << std::setw(6) <<  std::scientific << it->GroundStateWidth("neutron") << "\t"
                      << std::setw(6) <<  std::scientific << it->GroundStateWidth("proton") << "\t"
                      << std::setw(6) <<  std::scientific << it->GroundStateWidth("alpha") << "\t"
                      << std::setw(6) <<  std::scientific << it->TotalWidth("gamma") << "\t"
                      << std::setw(6) <<  std::scientific << it->TotalWidth("neutron") << "\t"
                      << std::setw(6) <<  std::scientific << it->TotalWidth("proton") << "\t"
                      << std::setw(6) <<  std::scientific << it->TotalWidth("alpha") << "\t"                       
                      << std::setw(6) <<  std::scientific << it->TotalWidth("all") << "\t"     
                      << std::fixed << std::endl;
        }
    }

    void WriteCompoundStates(CompoundStates * compound, std::string file){

    }


    std::string pTypeStringFromString(std::string &reactionString){
        std::string massNumberString=massNumberStringFromString(reactionString);
        std::string atomicNumberString=atomicNumberStringFromString(reactionString);

        reactionString.erase(0,massNumberString.length());
        reactionString.erase(0,atomicNumberString.length());
      
        std::string projectileString;
        for(unsigned int i = 0; i<reactionString.length(); i++) {
            std::string nextChar(reactionString,i,1);
            if(nextChar=="+") continue;
            else projectileString+=nextChar;
        }

        return projectileString;
    }

    int pTypeIntFromString(std::string &reactionString){
        std::string projectileString = pTypeStringFromString(reactionString);

        if(projectileString=="g")
            return 0;
        else if(projectileString=="n") 
            return 1;      
        else if(projectileString=="p")
            return 2;
        else if(projectileString=="a")
            return 3;
        
        /** Return "neutron" by default*/
        return 1;
    }

    bool fexists(const char *filename) {
        std::ifstream ifile(filename);
        return (bool)ifile;
    }

    std::string massNumberStringFromString(std::string &isotopeString){
        std::string massNumberString;

        for(unsigned int i = 0; i<isotopeString.length(); i++) {
      	    std::string nextChar(isotopeString,i,1);
            std::istringstream stm(nextChar);
            int nextDigit;
            if(!(stm>>nextDigit)) break;
            else massNumberString+=nextChar;
      }
      return massNumberString;
    }

    int massNumberIntFromString(std::string &isotopeString){
        std::string massNumberString = massNumberStringFromString(isotopeString);
        if(massNumberString.length()>0)
            return atoi(massNumberString.c_str());
        else
            return 0;
    }

    std::string atomicNumberStringFromString(std::string &reactionString){
        std::string massNumberS = massNumberStringFromString(reactionString);
        reactionString.erase(0,massNumberS.length());
        std::string atomicNumberString;
        for(unsigned int i = 0; i<reactionString.length(); i++) {
            std::string nextChar(reactionString,i,1);
            if(nextChar=="+") break;
                else atomicNumberString+=nextChar;
            }        
        return atomicNumberString;
    }

    int atomicNumberIntFromString(std::string &reactionString){
        std::string atomicNumberString = atomicNumberStringFromString(reactionString);
        if(NuclearMass::FindZ(atomicNumberString) != -1) 
            return NuclearMass::FindZ(atomicNumberString);
        else
            return 0;
    }

    void PrintHelp(){
        std::cout << std::endl << "\tSyntax:        sapphire compound <xxY+z> <energyFile>" << std::endl;        
        std::cout << std::endl << "\t               sapphire compound <xxY+z> <energyFile> <spinParity>" << std::endl;
	    std::cout << std::endl;        
    }

    void parseJPiString(double & J, int & Pi, std::string & jPiString){
        bool foundDelimiter = false;
        bool goodPi = false;
        bool goodJ = false;
        
        std::string firstString,secondString,parityString,delimiterString;

        for(int i = 0;i<jPiString.length();i++) {

        	std::string nextChar(jPiString,i,1);
            

        	if(nextChar == "/" || nextChar == ".") {
        	  foundDelimiter=true;
        	  delimiterString=nextChar;
        	  continue;

        	} else if(nextChar=="+" || nextChar=="-") {
        	    parityString = nextChar;
                break;
            }
    
            std::istringstream stm(nextChar);
            int digit;
            
            if(stm>>digit) {
                if(!foundDelimiter){
                    firstString+=nextChar;
                }
                else{
                    secondString+=nextChar;
                }
            }
        }

        if(parityString.length()>0) {
          if(parityString=="-") Pi=-1;
          else if(parityString=="+") Pi=1;
          goodPi=true;
        }

        if(firstString.length()>0) {
          	if(foundDelimiter&&delimiterString=="/") {
         
          	    if(secondString.length()>0){
                    J = atof(firstString.c_str())/atof(secondString.c_str());
         
                }
          	    else{
                    J = atof(firstString.c_str());
                }

          	} else if(foundDelimiter&&delimiterString==".") {
          	  firstString += '.'+secondString;
          	  J = atof(firstString.c_str());

          	} else {
                J = atof(firstString.c_str());
            }

          	double intPart;
            
            if(modf(J*2.,&intPart)==0.) goodJ=true; 

          }
    }


    double spinDoubleFromString(std::string &jPiString){
      
      double J;
      int Pi;

      parseJPiString(J,Pi, jPiString);
      
      return J;
    }

    int parityIntFromString(std::string &jPiString){
      
      double J;
      int Pi;

      parseJPiString(J,Pi, jPiString);
      
      return Pi;
    }
}