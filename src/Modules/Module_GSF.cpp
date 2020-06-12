#include "Modules/Module_GSF.h"
#include "NuclearMass.h"
#include "GammaStrength/GammaTransmissionFunc.h"

#include <iostream>
#include <string>
#include <sstream>



namespace Module_GSF{

    void Go(int argc, char *argv[]){
        if(argc < 6){
            std::cout << std::endl << "Not enough input parameters ..." << std::endl;
            PrintHelp();
            exit(1);
        }

        std::string isotope(argv[2]);

        GetGSF(isotope, atoi(argv[3]), atoi(argv[4]), atoi(argv[5]));

    }

    void GetGSF(std::string isotope,int e1, int m1, int e2){
        //std::cout << "Calling GetGS" << std::endl;
        ///std::cout << std::endl  << e1 << "\t"
        //                        << m1 << "\t"
        //                        << e2 << "\t"
        //                        << std::endl;
        
        int charge = atomicNumberIntFromString(isotope);
        int mass = massNumberIntFromString(isotope);

        GammaTransmissionFunc::SetEGDRType(e1);
        GammaTransmissionFunc::SetMGDRType(m1);
        GammaTransmissionFunc::SetEGQRType(e2);

        GammaTransmissionFunc* e1p;
        GammaTransmissionFunc* m1p;
        GammaTransmissionFunc* e2p;

        GammaTransmissionFunc* e1s  = GammaTransmissionFunc::CreateGammaTransmissionFunc(charge, mass, 0, 1, 1,-1, 0, 0, 0, 0, e1p, 0);
        GammaTransmissionFunc* m1s  = GammaTransmissionFunc::CreateGammaTransmissionFunc(charge, mass, 0, 1, 1, 1, 1, 0, 0, 0, m1p, 0);
        GammaTransmissionFunc* e2s  = GammaTransmissionFunc::CreateGammaTransmissionFunc(charge, mass, 0, 1, 2, 1, 2, 0, 0, 0, e2p, 0);
        
        std::cout.precision(2);
        std::cout << std::endl  << "Energy\t\t" << "f(E1)\t\t" << "f(M1)\t\t" << "f(E2)\t\t" << "T(E1)\t\t" << "T(M1)\t\t" << "T(E2)\t\t" << std::endl;

        for(double energy =0.5; energy <= 20.0; energy+=0.5){
            
            std::cout << std::scientific << energy << "\t"
                                    << e1s->CalcStrengthFunction(energy) << "\t"
                                    << m1s->CalcStrengthFunction(energy) << "\t"
                                    << e2s->CalcStrengthFunction(energy) << "\t"
                                    << e1s->operator()(energy) << "\t"
                                    << m1s->operator()(energy) << "\t"
                                    << e2s->operator()(energy) << "\t"
                                    << std::endl;
        }
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

    std::string atomicNumberStringFromString(std::string &isotopeString){
        std::string massNumberString = massNumberStringFromString(isotopeString);
        std::string atomicNumberString = isotopeString.substr(massNumberString.length());
        return atomicNumberString;
    }

    int atomicNumberIntFromString(std::string &isotopeString){
        std::string atomicNumberString = atomicNumberStringFromString(isotopeString);
        int Z = NuclearMass::FindZ(atomicNumberString);
        return Z;
    }

    void PrintHelp(){
        std::cout << std::endl << "Module: gsf" << std::endl; 
        std::cout << std::endl << "\tSyntax:        sapphire gsf E1Model M1Model E2Model" << std::endl;        
	    std::cout << std::endl;        
    }
}