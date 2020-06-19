#include "Modules/Module_Compound.h"
#include "Databases/NuclearMass.h"
#include "CrossSection.h"


#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <chrono>



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
        
        CalcWidths(reaction, atof(argv[3]));

        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << std::endl << "Total calculation time: " << elapsed_seconds.count() << "s\n";

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
            xs->CalcCompoundWidth();
        }
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
        std::cout << std::endl << "Module: compound" << std::endl; 
        std::cout << std::endl << "\tSyntax:        sapphire compound <xxY+z> <energyFile>" << std::endl;        
	    std::cout << std::endl;        
    }
}