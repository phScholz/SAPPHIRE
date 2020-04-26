/**
 * @file Module_RandomScheme.h
 * @date 2020-04-25
 * @brief Module to create random level schemes and calculate features.
 * 
 */

#include "RandomScheme.h"
#include "SapphireInput.h"
#include "NuclearLevels.h"

#include <iostream>
#include <chrono>
#include <vector>
#include <stdexcept>

namespace Module_RandomScheme{
    void PrintHelp(){
        std::cout  << "\tSyntax:        sapphire random <options>" << std::endl;        
	    std::cout << std::endl << "Options:" << std::endl;
        std::cout << std::endl;
        std::cout << "\tInputFile             " << "- determine input parameters from InputFile and run calculations." << std::endl << std::endl;
        std::cout << "\tcreate Z A Emax       " << "- Create a random level scheme for an isotope:" << std::endl;
        std::cout << "\t                      " << "   Z    = charge" << std::endl;
        std::cout << "\t                      " << "   A    = mass" << std::endl;
        std::cout << "\t                      " << "   Emax = maximum excitation energy" << std::endl << std::endl;
        std::cout << "\textend Z A Emax       " << "- Extend a level scheme for an isotope with random levels:" << std::endl;
        std::cout << "\t                      " << "   Z    = charge" << std::endl;
        std::cout << "\t                      " << "   A    = mass" << std::endl;
        std::cout << "\t                      " << "   Emax = maximum excitation energy" << std::endl;
        std::cout << std::endl; 
    }

    void Create(int argc, char* argv[]){
        if(argc < 6){
            throw std::invalid_argument("Invalid Input Parameters.");
        }
        int Z=atoi(argv[3]);
        int A=atoi(argv[4]);
        double E=atof(argv[5]);

        RandomScheme* scheme = new RandomScheme();

        
        std::cout << "Creating Level Scheme for" <<std::endl;
        std::cout << "\tZ = "<< Z << std::endl;
        std::cout << "\tA = "<< A << std::endl;
        std::cout << "\tE = "<< E << std::endl << std::endl;
        
        try
        {
            scheme->CreateRandomScheme(Z,A,0.0,E);
        }
        catch(const std::exception &e)
        {
            std::cerr << e.what() << '\n';
        }
        
        std::cout << "Printing Level Scheme ..." <<std::endl;
        try
        {
            scheme->PrintRandomScheme(E);    
        }
        catch(const std::exception& e)
        {
            std::cerr << e.what() << '\n';
        }

        delete scheme;       
    }

    void Extend(int argc, char* argv[]){
        if(argc < 6){
            throw std::invalid_argument("Invalid Input Parameters.");
        }
        int Z=atoi(argv[3]);
        int A=atoi(argv[4]);
        double E=atof(argv[5]);

        RandomScheme* scheme = new RandomScheme();

        
        std::cout << "Extending Level Scheme for" <<std::endl;
        std::cout << "\tZ = "<< Z << std::endl;
        std::cout << "\tA = "<< A << std::endl;
        std::cout << "\tE = "<< E << std::endl << std::endl;
        
        try
        {
            scheme->ExtendRandomScheme(Z,A,E);
        }
        catch(const std::exception &e)
        {
            std::cerr << e.what() << '\n';
        }
        
        std::cout << "Printing Level Scheme ..." <<std::endl;
        try
        {
            scheme->PrintRandomScheme(E);    
        }
        catch(const std::exception& e)
        {
            std::cerr << e.what() << '\n';
        }

        delete scheme;       
    }

    void Go(int argc,char* argv[]){
        /**
        *   A clock is started at the beginning of Go() 
        *   to measure the total calculation time. 
        */
        auto start = std::chrono::steady_clock::now();
        std::cout << std::endl << "Module: random" << std::endl;
        std::cout << std::endl;

        if(argc < 3){PrintHelp(); exit(0);}

        std::string mode(argv[2]);

        if(mode=="create")
        {
            try
            {
                Create(argc,argv);
            }
            catch(std::invalid_argument &e)
            {
                PrintHelp();
            }
        }
        else if (mode=="extend")
        {
            try
            {
                Extend(argc,argv);
            }
            catch(std::invalid_argument &e)
            {
                PrintHelp();
            }
        }
        else
        {
            PrintHelp();
        }
        
        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << std::endl << "Total calculation time: " << elapsed_seconds.count() << "s\n";
    }


}