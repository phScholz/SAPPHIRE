#include "Modules/Module_NLD.h"
#include "Databases/NuclearMass.h"
#include "LevelDensity/LevelDensity.h"
#include "LevelDensity/RauscherLevelDensity.h"
#include "LevelDensity/LevelDensityHFB_BSk14.h"

#include <iostream>
#include <string>
#include <sstream>



namespace Module_NLD{

    void Go(int argc, char *argv[]){
        if(argc < 4){
            std::cout << std::endl << "Not enough input parameters ..." << std::endl;
            PrintHelp();
            exit(1);
        }

        std::string isotope(argv[2]);
        
        std::cout << std::endl << "Nuclear Level Density for " << isotope << std::endl;

        GetNLD(isotope, atoi(argv[3]));

    }

    void GetNLD(std::string isotope, int nldmodel){
        int charge = atomicNumberIntFromString(isotope);
        int mass = massNumberIntFromString(isotope);

        switch(nldmodel){
            case 0:
                std::cout << "NLD model:\t" << "RauscherLevelDensity" << std::endl;
                Print<RauscherLevelDensity>(isotope);
                break;
            case 1:
                std::cout << "NLD model:\t" << "HFB_BSk14" << std::endl;
                Print<LevelDensityHFB_BSk14>(isotope);
                break;
            default:
                std::cout << "!!!Wrong input for NLD model !!!" << std::endl;
                exit(1);
        }
    }

    template <class T>
    void Print(std::string isotope){
        int charge = atomicNumberIntFromString(isotope);
        int mass = massNumberIntFromString(isotope);

        if(mass%2==0){
            LevelDensity * p0 = new T(charge, mass, 0, 1 );
            LevelDensity * n0 = new T(charge, mass, 0, -1 );
            LevelDensity * p1 = new T(charge, mass, 1, 1 );
            LevelDensity * n1 = new T(charge, mass, 1, -1 );
            LevelDensity * p2 = new T(charge, mass, 2, 1 );
            LevelDensity * n2 = new T(charge, mass, 2, -1 );
            LevelDensity * p3 = new T(charge, mass, 3, 1 );
            LevelDensity * n3 = new T(charge, mass, 3, -1 );

            std::cout << std::endl  << "Energy\t "<< "p(0+) \t"
                                    << "p(0-) \t"
                                    << "p(1+) \t"
                                    << "p(1-) \t"
                                    << "p(2+) \t"
                                    << "p(2-) \t"
                                    << "p(3+) \t"
                                    << "p(3-) \t"
                                    << std::endl;
            
            std::cout.precision(2);

            for(double energy=0.0; energy <=20.0; energy+=0.5){
                std::cout   << std::scientific << energy << " "
                            << p0->operator()(energy) << " " << n0->operator()(energy) << " "
                            << p1->operator()(energy) << " " << n1->operator()(energy) << " "
                            << p2->operator()(energy) << " " << n2->operator()(energy) << " "
                            << p3->operator()(energy) << " " << n3->operator()(energy) << " "
                            << std::endl;
            }
        }

        if(mass%2==1){
            LevelDensity * p0 = new T(charge, mass, 0.5, 1 );
            LevelDensity * n0 = new T(charge, mass, 0.5, -1 );
            LevelDensity * p1 = new T(charge, mass, 1.5, 1 );
            LevelDensity * n1 = new T(charge, mass, 1.5, -1 );
            LevelDensity * p2 = new T(charge, mass, 2.5, 1 );
            LevelDensity * n2 = new T(charge, mass, 2.5, -1 );
            LevelDensity * p3 = new T(charge, mass, 3.5, 1 );
            LevelDensity * n3 = new T(charge, mass, 3.5, -1 );

            std::cout << std::endl  << "Energy\t "<< "p(1/2+)\t"
                                    << "p(1/2-)\t"
                                    << "p(3/2+)\t"
                                    << "p(3/2-)\t"
                                    << "p(5/2+)\t"
                                    << "p(5/2-)\t"
                                    << "p(7/2+)\t"
                                    << "p(7/2-)\t"
                                    << std::endl;
            
            std::cout.precision(2);

            for(double energy=0.0; energy <=20.0; energy+=0.5){
                std::cout   << std::scientific << energy << " "
                            << p0->operator()(energy) << " " << n0->operator()(energy) << " "
                            << p1->operator()(energy) << " " << n1->operator()(energy) << " "
                            << p2->operator()(energy) << " " << n2->operator()(energy) << " "
                            << p3->operator()(energy) << " " << n3->operator()(energy) << " "
                            << std::endl;
            }
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
        std::cout << std::endl << "Module: nld" << std::endl; 
        std::cout << std::endl << "\tSyntax:        sapphire nld NLDModel" << std::endl;        
	    std::cout << std::endl;        
    }
}