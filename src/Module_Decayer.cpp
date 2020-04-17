#include "Module_Decayer.h"
#include <iostream>
#include <fstream>
#include "SapphireMPITypes.h"
#include "SapphireInput.h"
#include <string>
#include "NuclearMass.h"
#include <chrono>
#include <sstream>

namespace Module_Decayer{

    std::string massNumberStringFromString(std::string isotopeString){
        std::string massNumberString;

        for(int i = 0; i<isotopeString.length(); i++) {
      	    std::string nextChar(isotopeString,i,1);
            std::istringstream stm(nextChar);
            int nextDigit;
            if(!(stm>>nextDigit)) break;
            else massNumberString+=nextChar;
      }
      return massNumberString;
    }

    int massNumberIntFromString(std::string isotopeString){
        std::string massNumberString = massNumberStringFromString(isotopeString);
        if(massNumberString.length()>0)
            return atoi(massNumberString.c_str());
        else
            return 0;
    }

    std::string atomicNumberStringFromString(std::string isotopeString){
        std::string massNumberString = massNumberStringFromString(isotopeString);
        std::string atomicNumberString = isotopeString.substr(massNumberString.length());
        return atomicNumberString;
    }

    int atomicNumberIntFromString(std::string isotopeString){
        std::string atomicNumberString = atomicNumberStringFromString(isotopeString);
        int Z = NuclearMass::FindZ(atomicNumberString);
        return Z;
    }

    void Go(int argc,char *argv[]){
        auto start = std::chrono::steady_clock::now();
        std::cout << "Module: decayer" << std::endl;
        std::cout << std::endl;

        if(argc < 3){
            printHelp();
            exit(0);
        }
        
        if(fexists(argv[2])){              
            SapphireInput* Input = new SapphireInput();    
            std::string str(argv[2]);
            Input->printIntputFile(str);
            Input->ReadInputFile(str);
            Input->printIntputParameters();

            delete Input;
        }else{
            std::cout << "No valid input file was given." << std::endl;
            exit(1);
        }

        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << std::endl << "Total calculation time: " << elapsed_seconds.count() << "s\n";
    }
    
    void Run(int argc,char *argv[]){}

    bool fexists(const char *filename) {
        std::ifstream ifile(filename);
    return (bool)ifile;
    }

    void printHelp(){
        std::cout  << "\tSyntax:        sapphire decayer <options>" << std::endl;        
	    std::cout << std::endl << "Options:" << std::endl;
        std::cout << std::endl;
        /* std::cout << "\tAX+y           - reaction string, e.g. 60Fe+p, running calculations with default settings." << std::endl;
        std::cout << "\tInputFile      - determine input parameters from InputFile and run calculations." << std::endl;
        std::cout << std::endl; */
    }
}