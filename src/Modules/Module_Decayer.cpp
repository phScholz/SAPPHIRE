/**
 * @file Sapphire.cpp
 * @brief Entry point for the decayer module. * 
 */


#include "Modules/Module_Decayer.h"
#include "SapphireMPITypes.h"
#include "SapphireInput.h"
#include <string>
#include "Databases/NuclearMass.h"
#include <chrono>
#include "Decayer/Decayer.h"
#include "omp.h"
#include "Decayer/DecayController.h"
#include "Decayer/DecayResults.h"
#include "Progressbar.h"
#include "ParticleTransmissionFunc.h"
#include "GammaStrength/GammaTransmissionFunc.h"
#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"
#include <ctime>
#include "SPDistribution.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"

namespace Module_Decayer{

    bool ReadDist(std::string file, SPDistribution &dist){
        std::cout << std::endl << "Reading Spin-Parity Distribution ..." << std::endl;
        /**
         * 1. Check whether the SPDistribution object is empty. If not, return false.
         */
        if (dist.distribution == nullptr){

            /**
             * 2. Attempting to read distribution file. If not possible, return false.
             */
            std::ifstream in(file.c_str());
            if(!in){
                std::cout << std::endl << "Cannot read distribution file ... " << file << std::endl;
                return false;
            }

            /**
             * 3. Read content of distribution file line-by-line and store it into the SPDistribution object.
             */
            std::string line;
            SPPopulation dummy;
            std::vector<SPPopulation> *dummyVector = new std::vector<SPPopulation>;

            double spin;
            double pop;
            int parity;

            while(!in.eof()){
                std::getline(in,line);
                std::istringstream lineStream(line);
                lineStream >> spin >> parity >> pop;
                std::cout << "\t" << spin << "\t" << parity << "\t" << pop << std::endl;
                dummy.Spin(spin);
                dummy.Parity(parity);
                dummy.Pop(pop);                
                dummyVector->push_back(dummy);                
            }

            in.close();
            dist.distribution = dummyVector;
            dist.Normalize();
            dist.PrintPopulation();
            
        }
        else{
            std::cout << std::endl << "Distribution object was not initialized!!" << std::endl;
            return false;
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
            
            if (Input.DisFile().length()>0){
                std::cout << std::endl << "Running Decayer for distribution in ... " << Input.DisFile() << std::endl;
                RunDist(Input);
            }
            else{
                std::cout << std::endl << "Running Decayer for resonances with J = " << Input.Spin() << " and Pi = " << Input.Parity() << std::endl;
                RunSingle(Input);
            }
            
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

    void RunDist(const SapphireInput &input){
        /**
         * 1. Attempting to read the spin-parity distribution from file. If this fails, the code exits.
         */
        SPDistribution spinDist;
        if(!ReadDist(input.DisFile(), spinDist)){
            std::cout << std::endl << "Cannot load spin-parity distribution!!!" << std::endl;
            exit(1);
        }

        /**
        * 2. The parameters needed to be initialized for the old Sapphire code
        * are obtained from the SapphireInput object passed by reference to the Module_Decayer::Run()
        * method.
        */
        Decayer::SetCrossSection(false);
        input.SetInputDecayer();
        input.SetInputParticleTransmission();
        input.SetInputGammaTransmission();

        int chunkSize = input.ChunkSize();
        int events = input.Events();
        
        std::string isotopeString(input.Isotope());
        int A = massNumberIntFromString(isotopeString);
        int Z = atomicNumberIntFromString(isotopeString);
        int Pi = input.Parity();
        
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
        *   3. Because of the "segmentation fault on more than 10 threads" bug,
        *   the system is asked via omp_get_max_threads() what the maximum numbers of threads are.
        *   If this value is larger than 10, then the maximum number of threads used for the calculation
        *   is fixed at 10 to prevent the "segmentation fault" bug. 
        *   Once this issue is fixed, this can be removed.
        */
        if(omp_get_max_threads() > 10) omp_set_num_threads(10);

        /**
        * 4. If the number of events is not a multiple of the chunkSize,
        * then the number of remainder is calculated from the modulo.
        * The number of maximum chunks is then derived from the ratio (events-remainder)/chunkSize.
        */
        int remainder = events%chunkSize;
        int chunks = (events-remainder)/chunkSize;

        /**
        * 5. For the Monte-Carlo decay, we need to initialize a random generator.
        *    We are using here the boost::random library and initialize it with the current time.
        */

        std::cout << std::endl << "Starting Decay Simulation..." << std::endl;

        //Select random number generator
        typedef boost::mt19937 base_generator_type;
    
        //Use the system time as seed
        base_generator_type generator(time(0));
    
        // Define a uniform random number distribution which produces "double"
        // values between 0 and 1 (0 inclusive, 1 exclusive).
        boost::uniform_real<> uni_dist(0,1);
        boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(generator, uni_dist);

        // Initialize the Decay results object
        DecayResults* results = NULL;
        if(events>1) results = new DecayResults(Z,A,lowEnergy,highEnergy,suffixNo);
        
        
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
            double energy=0;

            std::vector<DecayController *> controllerVector(numInChunk, nullptr);
                      
            for(int k = 0; k<numInChunk; k++){
                /**
                *   We randomly draw the initial energy.
                */
                energy = (lowEnergy==highEnergy) ? lowEnergy : lowEnergy+(highEnergy-lowEnergy)*uni();

                /**
                *   We randomly draw the initial spin and parity
                */

                double J;
                int Pi;
                double rdmCDF = uni();

                for(std::vector<SPPopulation>::const_iterator it = spinDist.distribution->begin(); it!= spinDist.distribution->end(); ++it){
                    if(it->Cdf() >= rdmCDF){
                        J = it->Spin();
                        Pi = it->Parity();
                        break;
                    }
                }
                
                
                //std::cout << "\t" << J << "\t" << Pi << std::endl;

                if(preEq)
                {
	               controllerVector.at(k)= new DecayController(Z,A,J,Pi,energy,numNuParticles,numNuHoles,numPiParticles,numPiHoles);
                } 
                else
                {
                    controllerVector.at(k) = new DecayController(Z,A,J,Pi,energy);
                }
            }

            pg.start(numInChunk); 

            #pragma omp parallel for
            for(int j = 0;j<numInChunk;j++){
                
                pg.update(j);

                DecayController* controller = controllerVector.at(j);

                double neutronEntranceWidth = 0.;
                double protonEntranceWidth = 0.;
                double gammaEntranceWidth = 0.;
                double alphaEntranceWidth = 0.;
                double neutronTotalWidth = 0.;
                double protonTotalWidth = 0.;
                double gammaTotalWidth = 0.;
                double alphaTotalWidth = 0.;

                controller->Decay(neutronEntranceWidth,protonEntranceWidth,alphaEntranceWidth,gammaEntranceWidth,neutronTotalWidth,protonTotalWidth,alphaTotalWidth,gammaTotalWidth); 
                
                #pragma omp critical
                    chunkResults[j] = std::pair<DecayData,std::vector<DecayProduct>>(DecayData(controller->Energy(), controller->Spin(), controller->Parity(), neutronEntranceWidth,protonEntranceWidth, alphaEntranceWidth,gammaEntranceWidth, neutronTotalWidth,protonTotalWidth, alphaTotalWidth,gammaTotalWidth),controller->DecayProducts());
                
                #pragma omp critical
                    if(events<=10) controller->PrintDecays();
                //delete controller;
            }

            pg.update(numInChunk);

            if(events>1){
                std::cout << std::endl << "Writing ROOT Tree..." << std::endl;                
                results->AddResults(chunkResults);
            }

            while(!controllerVector.empty()) delete controllerVector.back(), controllerVector.pop_back();            
        }

        if(events>1)
        {
            results->WriteNCloseFile();
            delete results;
        }
    }
    
    void RunSingle(const SapphireInput &input){
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
        int events = input.Events();
        
        std::string isotopeString(input.Isotope());
        int A = massNumberIntFromString(isotopeString);
        int Z = atomicNumberIntFromString(isotopeString);
        int Pi = input.Parity();
        
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
            double energy=0;

            std::vector<DecayController *> controllerVector(numInChunk, nullptr);

            //Select random number generator
            typedef boost::mt19937 base_generator_type;

            //Use the system time as seed
            base_generator_type generator(time(0));

            // Define a uniform random number distribution which produces "double"
            // values between 0 and 1 (0 inclusive, 1 exclusive).
            boost::uniform_real<> uni_dist(0,1);

            boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(generator, uni_dist);

            for(int k = 0; k<numInChunk; k++){
                energy = (lowEnergy==highEnergy) ? lowEnergy : lowEnergy+(highEnergy-lowEnergy)*uni();
                std::cout << std::endl << energy << std::endl;

                if(preEq)
                {
	               controllerVector.at(k)= new DecayController(Z,A,J,Pi,energy,numNuParticles,numNuHoles,numPiParticles,numPiHoles);
                } 
                else
                {
                    controllerVector.at(k) = new DecayController(Z,A,J,Pi,energy);
                }
            }

            pg.start(numInChunk); 

            #pragma omp parallel for
            for(int j = 0;j<numInChunk;j++){
                
                pg.update(j);

                DecayController* controller = controllerVector.at(j);

                double neutronEntranceWidth = 0.;
                double protonEntranceWidth = 0.;
                double gammaEntranceWidth = 0.;
                double alphaEntranceWidth = 0.;
                double neutronTotalWidth = 0.;
                double protonTotalWidth = 0.;
                double gammaTotalWidth = 0.;
                double alphaTotalWidth = 0.;

                controller->Decay(neutronEntranceWidth,protonEntranceWidth,alphaEntranceWidth,gammaEntranceWidth,neutronTotalWidth,protonTotalWidth,alphaTotalWidth,gammaTotalWidth); 
                
                #pragma omp critical
                    chunkResults[j] = std::pair<DecayData,std::vector<DecayProduct>>(DecayData(controller->Energy(), J, Pi, neutronEntranceWidth,protonEntranceWidth, alphaEntranceWidth,gammaEntranceWidth, neutronTotalWidth,protonTotalWidth, alphaTotalWidth,gammaTotalWidth),controller->DecayProducts());
                
                #pragma omp critical
                    if(events<=10) controller->PrintDecays();
                //delete controller;
            }

            pg.update(numInChunk);

            if(events>1){
                std::cout << std::endl << "Writing ROOT Tree..." << std::endl;                
                results->AddResults(chunkResults);
            }

            while(!controllerVector.empty()) delete controllerVector.back(), controllerVector.pop_back();            
        }

        if(events>1)
        {
            results->WriteNCloseFile();
            delete results;
        }
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