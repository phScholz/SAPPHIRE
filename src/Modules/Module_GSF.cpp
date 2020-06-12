#include "Modules/Module_GSF.h"
#include "Databases/NuclearMass.h"
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
        
        std::cout << std::endl << "Gamma-Strength Functions for " << isotope << std::endl;

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

        switch (e1)
        {
        case 0:
            std::cout << "E1 Strength Function:\t" << "Brink-AxelGSF" << std::endl;
            break;

        case 1:
            std::cout << "E1 Strength Function:\t" << "KopeckyUhlGSF" << std::endl;
            break;
        
        case 2:
            std::cout << "E1 Strength Function:\t" << "McCullaghGSF" << std::endl;
            break;
        
        case 3:
            std::cout << "E1 Strength Function:\t" << "D!MQRPA" << std::endl;
            break;
        
        default:
            std::cout << "!!!Wrong Input for E1 Strength Function!!!" << std::endl;
            exit(1);
            break;
        }

        switch (m1)
        {
        case 0:
            std::cout << "M1 Strength Function:\t" << "Brink-AxelGSF" << std::endl;
            break;

        //case 1:
        //    std::cout << "E1 Strength Function:" << "KopeckyUhlGSF" << std::endl;
        //    break;
        
        case 2:
            std::cout << "M1 Strength Function:\t" << "McCullaghGSF" << std::endl;
            break;
        
        case 3:
            std::cout << "M1 Strength Function:\t" << "D!MQRPA" << std::endl;
            break;
        
        default:
            std::cout << "!!!Wrong Input for M1 Strength Function!!!" << std::endl;
            exit(1);
            break;
        }

        switch (e2)
        {
        case 0:
            std::cout << "E2 Strength Function:\t" << "Brink-AxelGSF" << std::endl;
            break;

        //case 1:
        //    std::cout << "E1 Strength Function:" << "KopeckyUhlGSF" << std::endl;
        //    break;
        
        case 2:
            std::cout << "E2 Strength Function:\t" << "McCullaghGSF" << std::endl;
            break;
        
        //case 3:
        //    std::cout << "E1 Strength Function:" << "D!MQRPA" << std::endl;
        //    break;
        
        default:
            std::cout << "!!!Wrong Input for E2 Strength Function!!!" << std::endl;
            exit(1);
            break;
        }

        GammaTransmissionFunc* e1p;
        GammaTransmissionFunc* m1p;
        GammaTransmissionFunc* e2p;

        GammaTransmissionFunc* e1s  = GammaTransmissionFunc::CreateGammaTransmissionFunc(charge, mass, 0, 1, 1,-1, 0, 0, 0, 0, e1p, 0);
        GammaTransmissionFunc* m1s  = GammaTransmissionFunc::CreateGammaTransmissionFunc(charge, mass, 0, 1, 1, 1, 1, 0, 0, 0, m1p, 0);
        GammaTransmissionFunc* e2s  = GammaTransmissionFunc::CreateGammaTransmissionFunc(charge, mass, 0, 1, 2, 1, 2, 0, 0, 0, e2p, 0);
        
        std::cout.precision(2);
        std::cout << std::endl  << "Energy\t\t" << "f(E1)\t\t" << "f(M1)\t\t" << "f(E2)\t\t" << "T(E1)\t\t" << "T(M1)\t\t" << "T(E2)\t\t" << std::endl;

        for(double energy =0.0001; energy <= 25.0; energy+=0.5){
            
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
        std::cout << std::endl << "\tSyntax:        sapphire gsf Isotope E1Model M1Model E2Model" << std::endl;        
	    std::cout << std::endl;        
    }
}