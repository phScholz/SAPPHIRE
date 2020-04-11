#include "Module_CrossSection.h"
#include <iostream>
#include <fstream>
#include <sstream>

#include "NuclearMass.h"

namespace Module_CrossSection{
    void Go(int argc,char *argv[]){
       printHelp();
    }

    void Run(int argc,char *argv[]){}

    bool fexists(const char *filename) {
        std::ifstream ifile(filename);
    return (bool)ifile;
    }

    int pTypeFromString(const char* projectileString){
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