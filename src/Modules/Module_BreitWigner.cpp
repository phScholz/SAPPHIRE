#include "Modules/Module_BreitWigner.h"
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




namespace Module_BreitWigner{
    void Go(int argc, char *argv[]){
        auto start = std::chrono::steady_clock::now();
        std::cout << std::endl << "Module: breitWigner" << std::endl;
        std::cout << std::endl;
        if(argc < 3){
            std::cout << std::endl << "Not enough input parameters ..." << std::endl;
            PrintHelp();
            exit(1);
        }

        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << std::endl << "Total calculation time: " << elapsed_seconds.count() << "s\n";


    }

    void PrintHelp(){
        std::cout << std::endl << "\tSyntax:        sapphire breitWigner <xxY+z> <partialWidthsFile> <energyFile>" << std::endl;        
        std::cout << std::endl << "\t               sapphire breitWigner <InputFile>" << std::endl;
	    std::cout << std::endl;
    }
}