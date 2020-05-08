/**
 * @file Sapphire.cpp
 * @brief Entry point for the decayer module. * 
 */


#include "Module_Decayer.h"
#include <iostream>
#include <fstream>
#include "SapphireMPITypes.h"
#include "SapphireInput.h"
#include <string>
#include "NuclearMass.h"
#include <chrono>
#include <sstream>
#include "Decayer.h"
#include "omp.h"
#include "DecayController.h"
#include "DecayResults.h"
#include "Progressbar.h"
#include "ParticleTransmissionFunc.h"
#include "GammaStrength/GammaTransmissionFunc.h"

namespace Module_Decayer{

    std::string massNumberStringFromString(std::string isotopeString){
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
        /**
        *   A clock is started at the beginning of Go() 
        *   to measure the total calculation time. 
        */
        auto start = std::chrono::steady_clock::now();
        std::cout << "Module: decayer" << std::endl;
        std::cout << std::endl;
        /** 
        * If too few arguments are given for the decayer, being "program module inputFile",
        * then the printHelp() method is called and the program exits.
        */
        if(argc < 3){
            printHelp();
            exit(0);
        }
        
        /**
        * In a next step it'll checked if the filename points to an actual file.
        * If yes, the calculation can start. If not, the program exits.
        */
        if(fexists(argv[2])){     
            /**
            * If the inputFile exists:
            * - an SapphireInput object is created
            * - the content of the inputFile is printed to std::cout via printInputFile()
            * - the inputFile is read by ReadInputFile()
            * - the parameters in the Sapphire Input Object are printed via printInputParameters()
            * - the Run() method is called, passing the SapphireInput object by reference
            */
            SapphireInput Input;    
            std::string str(argv[2]);
            Input.PrintIntputFile(str);
            Input.ReadInputFile(str);
            Input.PrintIntputParameters("Decayer");
            Run(Input);
            
        }else{
            std::cout << "No valid input file was given." << std::endl;
            exit(1);
        }

        /** 
        * At the end of the Go() method the clock is stopped and the total calculation time
        * is printed. 
        */
        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << std::endl << "Total calculation time: " << elapsed_seconds.count() << "s\n";
    }
    
    void Run(const SapphireInput &input){
        /**
        * 1. At first the parameters needed to be initialized for the old Sapphire code
        * are obtained from the SapphireInput object passed by reference to the Module_Decayer::Run()
        * method.
        */
        Decayer::SetCrossSection(false);
        input.SetInputDecayer();
        input.SetInputParticleTransmission();
        input.SetInputGammaTransmission();

        int chunkSize = input.ChunkSize();
        int A = massNumberIntFromString(input.Isotope());
        int Z = atomicNumberIntFromString(input.Isotope());
        int Pi = input.Parity();
        int events = input.Events();
        double J = input.Spin();
        double lowEnergy = input.LowEnergy();
        double highEnergy = input.HighEnergy();
        int suffixNo = input.Suffix();
        bool preEq = input.PreEq();
        int numPiHoles = 0;
        int numPiParticles =0;
        int numNuHoles =0;
        int numNuParticles =0;

        /**
        *   2. Because of the "segmentation fault on more than 10 threads" bug,
        *   the system is asked via omp_get_max_threads() what the maximum numbers of threads are.
        *   If this value is larger than 10, then the maximum number of threads used for the calculation
        *   is fixed at 10 to prevent the "segmentation fault" bug. 
        *   Once this issue is fixed, this can be removed.
        */
        if(omp_get_max_threads() > 10) omp_set_num_threads(10);

        /**
        *   3. In a next step, the input parameters for the Decayer are printed to std::ccout.
        */

        std::cout << "Input Values For Parent Nucleus:" << std::endl
  	    << std::setw(15) << "Z:"             << std::setw(12) << Z          
  	    << std::setw(0) << std::endl
  	    << std::setw(15) << "A:"             << std::setw(12) << A          
  	    << std::setw(0) << std::endl
  	    << std::setw(15) << "J:"             << std::setw(12) << J          
  	    << std::setw(0) << std::endl
  	    << std::setw(15) << "Pi:"            << std::setw(12) << Pi         
  	    << std::setw(0) << std::endl
  	    << std::setw(15) << "energy (low):"  << std::setw(12) << lowEnergy  
  	    << std::setw(0) << std::endl
  	    << std::setw(15) << "energy (high):" << std::setw(12) << highEnergy 
  	    << std::setw(0) << std::endl
  	    << std::setw(15) << "events:"        << std::setw(12) << events     
  	    << std::setw(0) << std::endl
        << std::setw(15) << "chunk size:"    << std::setw(12) << chunkSize     
  	    << std::setw(0) << std::endl << std::endl;

        /**
        * 4. If the number of events is not a multiple of the chunkSize,
        * then the number of remainder is calculated from the modulo.
        * The number of maximum chunks is then derived from the ratio (events-remainder)/chunkSize.
        */
        int remainder = events%chunkSize;
        int chunks = (events-remainder)/chunkSize;

        /**
        * 5. For the Monte-Carlo decay, an unsigned int randomSeed[12] array is initialized.
        */
        unsigned int randomSeed[12];

        std::cout << "Starting Decay Simulation..." << std::endl;

        DecayResults* results = NULL;
        if(events>1) results = new DecayResults(Z,A,J,Pi,lowEnergy,highEnergy,suffixNo);
        
    for(int i = 0;i<=chunks;i++) {

          
            int numInChunk = (i==chunks) ? remainder : chunkSize;
            if(numInChunk==0) continue;
            if(events>=chunkSize)
              std::cout << "Decay chunk " << i+1 << " of " << chunks << " started ..." << std::endl;
            else
              std::cout << "Decay of " << numInChunk << " nuclei started ..." << std::endl;
            std::vector<std::pair<DecayData,std::vector<DecayProduct> > > chunkResults;
            chunkResults.resize(numInChunk);
            
            /**
            *   Now a ProgressBar object is generated from Progressbar.h, initialized with numInChunk, which ist the 
            *   number of decays in the current chunk.
            */
            ProgressBar pg;
            pg.start(numInChunk); 

            
            #pragma omp parallel for
            for(int j = 0;j<numInChunk;j++) {
                
                pg.update(j);
                double energy = (lowEnergy==highEnergy) ? lowEnergy :
    
                lowEnergy+(highEnergy-lowEnergy)*double(rand_r(&randomSeed[omp_get_thread_num()]))/double(RAND_MAX);

                DecayController* controller;

                if(preEq) {
	               controller= new DecayController(Z,A,J,Pi,energy,numNuParticles,numNuHoles,numPiParticles,numPiHoles);
                } else controller = new DecayController(Z,A,J,Pi,energy);

                double neutronEntranceWidth = 0.;
                double protonEntranceWidth = 0.;
                double gammaEntranceWidth = 0.;
                double alphaEntranceWidth = 0.;
                double neutronTotalWidth = 0.;
                double protonTotalWidth = 0.;
                double gammaTotalWidth = 0.;
                double alphaTotalWidth = 0.;

                controller->Decay(neutronEntranceWidth,protonEntranceWidth,alphaEntranceWidth,gammaEntranceWidth,neutronTotalWidth,protonTotalWidth,alphaTotalWidth,gammaTotalWidth); 
                
                chunkResults[j] = std::pair<DecayData,std::vector<DecayProduct>>(DecayData(energy,neutronEntranceWidth,protonEntranceWidth, alphaEntranceWidth,gammaEntranceWidth, neutronTotalWidth,protonTotalWidth, alphaTotalWidth,gammaTotalWidth),controller->DecayProducts());

                if(events==1) controller->PrintDecays();
                delete controller;
            }
            
                        
            if(events>1){
                std::cout << std::endl << "Writing ROOT Tree..." << std::endl;                
                results->AddResults(chunkResults);
            }
            
        }       

        if(results) delete results;

    }

    bool fexists(const char *filename) {
        std::ifstream ifile(filename);
    return (bool)ifile;
    }

    void printHelp(){
        std::cout  << "\tSyntax:        sapphire decayer <options>" << std::endl;        
	    std::cout << std::endl << "Options:" << std::endl;
        std::cout << std::endl;
        std::cout << "\tInputFile      - determine input parameters from InputFile and run calculations." << std::endl;        std::cout << std::endl; 
    }
}