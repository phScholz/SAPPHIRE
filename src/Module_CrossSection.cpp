#include "Module_CrossSection.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "NuclearMass.h"
#include "CrossSection.h"

#include "Setup.cpp"

extern void Initialize();

namespace Module_CrossSection{

    typedef struct EntrancePairs {
        EntrancePairs(int Z,int A,int pType) {
            Z_=Z;
            A_=A;
            pType_=pType;
        };
        int Z_;
        int A_;
        int pType_;
    } EntrancePairs;

    void Go(int argc,char *argv[]){
        if(argc < 3){
            printHelp();
        }
        
        if(fexists(argv[2])){
            std::cout << "Now I would handle the input-File." << std::endl;
        }
        else{
            std::vector<EntrancePairs> entrancePairs;
            int A = massNumberIntFromString(argv[2]);
            int Z = atomicNumberIntFromString(argv[2]);
            int pType = pTypeIntFromString(argv[2]);
            std::string energyFile;
            bool forRates = false;
            int entranceState = 0;
            std::vector<int> exitStates(4,-1);
            


            CrossSection* xs = new CrossSection(Z,A,pType,energyFile,forRates,entranceState,exitStates);
            if(xs->IsValid())
            {
                xs->Calculate();
                xs->PrintCrossSections();
            }
            else
            {
                std::cout << "Could not calculate cross section." << std::endl;    
            }            
            delete xs;
        }
    }

    void Run(int argc,char *argv[]){}

    bool fexists(const char *filename) {
        std::ifstream ifile(filename);
    return (bool)ifile;
    }

    std::string pTypeStringFromString(std::string reactionString){
        massNumberString=massNumberStringFromString(reactionString);
        atomicNumberString=atomicNumberStringFromString(reactionString);

        reactionString.erase(0,massNumberString.length());
        reactionString.erase(0,atomicNumberString.length());
      
        std::string projectileString;
        for(int i = 0; i<reactionString.length(); i++) {
            std::string nextChar(reactionString,i,1);
            if(nextChar=="+") continue;
            else projectileString+=nextChar;
        }

        return projectileString;
    }

    int pTypeIntFromString(std::string reactionString){
        std::string projectileString = pTypeStringFromString(reactionString);

        if(projectileString=="g")
            return 0;
        else if(projectileString=="n") 
            return 1;      
        else if(projectileString=="p")
            return 2;
        else if(projectileString=="a")
            return 3;
    }

    std::string massNumberStringFromString(std::string reactionString){
        std::string massNumberS;
        for(int i = 0; i<reactionString.length(); i++) {
            std::string nextChar(reactionString,i,1);
            std::istringstream stm(nextChar);
            int nextDigit;
            if(!(stm>>nextDigit)) break;
                else massNumberS+=nextChar;
        }

        return massNumberS;
    }

    int massNumberIntFromString(std::string reactionString){
        std::string massNumberS = massNumberStringFromString(reactionString);
        if(massNumberS.length()>0)
            return atoi(massNumberS.c_str());
        else
            return 0;
    }

    std::string atomicNumberStringFromString(std::string reactionString){
        std::string massNumberS = massNumberStringFromString(reactionString);
        reactionString.erase(0,massNumberS.length());
        std::string atomicNumberString;
        for(int i = 0; i<reactionString.length(); i++) {
            std::string nextChar(reactionString,i,1);
            if(nextChar=="+") break;
                else atomicNumberString+=nextChar;
            }        
        return atomicNumberString;
    }

    int atomicNumberIntFromString(std::string reactionString){
        std::string atomin = atomicNumberStringFromString(reactionString);
        if(NuclearMass::FindZ(atomicNumberString) != -1) 
            return NuclearMass::FindZ(atomicNumberString);
        else
            return 0;
    }

    void printHelp(){
        std::cout << "sapphire reaction <options>" << std::endl;
        std::cout << std::endl;
        std::cout << "Options:" << std::endl;
        std::cout << std::endl;
        std::cout << " AX+y           - reaction string, e.g. 60Fe+p (not yet implemented)" << std::endl;
        std::cout << " InputFile      - InputFile (not yet implemented)" << std::endl;
    }
}