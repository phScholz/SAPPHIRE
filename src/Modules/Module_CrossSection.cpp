/**
 * @file Sapphire.cpp
 * @brief Entry point for the reaction module.
 * 
 * 
 * 
 */
#include "Modules/Module_CrossSection.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>

/*Time measurement*/
#include <chrono>
//#include <ctime>

#include "Databases/NuclearMass.h"
#include "CrossSection.h"
#include "SapphireInput.h"
#include "Decayer/Decayer.h"
#include "TransitionRateFunc.h"
#include "ParticleTransmissionFunc.h"
#include "GammaStrength/GammaTransmissionFunc.h"
#include "LevelDensity/LevelDensityTable.h"
#include "omp.h"
#include "Progressbar.h"


namespace Module_CrossSection{
 

    bool fexists(const char *filename) {
        std::ifstream ifile(filename);
        return (bool)ifile;
    }

    std::string pTypeStringFromString(std::string &reactionString){
        std::string massNumberString=massNumberStringFromString(reactionString);
        std::string atomicNumberString=atomicNumberStringFromString(reactionString);

        reactionString.erase(0,massNumberString.length());
        reactionString.erase(0,atomicNumberString.length());
      
        std::string projectileString;
        for(unsigned int i = 0; i<reactionString.length(); i++) {
            std::string nextChar(reactionString,i,1);
            if(nextChar=="+") continue;
            else projectileString+=nextChar;
        }

        return projectileString;
    }

    int pTypeIntFromString(std::string &reactionString){
        std::string projectileString = pTypeStringFromString(reactionString);

        if(projectileString=="g")
            return 0;
        else if(projectileString=="n") 
            return 1;      
        else if(projectileString=="p")
            return 2;
        else if(projectileString=="a")
            return 3;
        
        /** Return "neutron" by default*/
        return 1;
    }

    std::string massNumberStringFromString(std::string &reactionString){
        std::string massNumberS;
        for(unsigned int i = 0; i<reactionString.length(); i++) {
            std::string nextChar(reactionString,i,1);
            std::istringstream stm(nextChar);
            int nextDigit;
            if(!(stm>>nextDigit)) break;
                else massNumberS+=nextChar;
        }

        return massNumberS;
    }

    int massNumberIntFromString(std::string &reactionString){
        std::string massNumberS = massNumberStringFromString(reactionString);
        if(massNumberS.length()>0)
            return atoi(massNumberS.c_str());
        else
            return 0;
    }

    std::string atomicNumberStringFromString(std::string &reactionString){
        std::string massNumberS = massNumberStringFromString(reactionString);
        reactionString.erase(0,massNumberS.length());
        std::string atomicNumberString;
        for(unsigned int i = 0; i<reactionString.length(); i++) {
            std::string nextChar(reactionString,i,1);
            if(nextChar=="+") break;
                else atomicNumberString+=nextChar;
            }        
        return atomicNumberString;
    }

    int atomicNumberIntFromString(std::string &reactionString){
        std::string atomicNumberString = atomicNumberStringFromString(reactionString);
        if(NuclearMass::FindZ(atomicNumberString) != -1) 
            return NuclearMass::FindZ(atomicNumberString);
        else
            return 0;
    }

    void printHelp(){
        std::cout  << "\tSyntax:        sapphire reaction <options>" << std::endl;        
	    std::cout << std::endl << "Options:" << std::endl;
        std::cout << std::endl;
        std::cout << "\tAX+y           - reaction string, e.g. 60Fe+p, running with default settings." << std::endl;
        std::cout << "\tInputFile      - determine input parameters from InputFile and run calculations." << std::endl;
        std::cout << std::endl;
    }

    void readEntrancePairs(std::vector<EntrancePairs> * entrancePairs, std::string reactionFile){
        std::ifstream in(reactionFile.c_str());
        
        if(!in) {
	        std::cout << "Could not open " << reactionFile << " for reading." << std::endl;
	        exit(1);
        } else {
	        std::cout << "Reading nuclei from " << reactionFile << ".::" << std::endl << std::endl;

	    while(!in.eof()) {
            int Z,A,pType =0;
	        std::string line;
	        std::getline(in,line);
	            
            if(!in.eof()) {
	            std::istringstream stm(line);
	            if(stm >> Z >> A >> pType)
                {
                    entrancePairs->push_back(EntrancePairs(Z,A,pType));
                }
                else{
                    std::cout << "Cannot interpret line: " << line << std::endl;
                }
	        }
	    }
	                
        in.close();
        }
    }

    void PrintEntrancePairs(std::vector<EntrancePairs> & entrancePairs){
        for(auto it = std::begin(entrancePairs); it != std::end(entrancePairs); ++it){
            std::cout << "Z: " << it->Z_ << " A: " << it->A_ << " Particle Type: " << it->pType_ << std::endl;
        }
    }

    void Run(const SapphireInput & input){
        //Copy the input parameters to respective classes
        //I know, this is still not really intuitive ... but I am working on it.
        input.SetInputCrossSection();
        input.SetInputDecayer();
        input.SetInputTransitionRate();
        input.SetInputParticleTransmission();
        input.SetInputGammaTransmission();
        input.SetInputLevelDensity();

        
        std::vector<int> exitStates(4,-1);
        exitStates[0]=input.g_ExitStates();
        exitStates[1]=input.n_ExitStates();
        exitStates[2]=input.p_ExitStates();
        exitStates[3]=input.a_ExitStates();        

        

        
        std::vector<EntrancePairs> entrancePairs;
        readEntrancePairs(&entrancePairs,input.ReactionFile());
        std::vector<EntrancePairs>::iterator it;

        //PrintEntrancePairs(entrancePairs);

        
        for(it = std::begin(entrancePairs); it != std::end(entrancePairs); ++it){
            CrossSection* xs = new CrossSection(it->Z_,it->A_, it->pType_,input.EnergyFile(),input.CalcRates());
            if(xs->IsValid())
            {
                if(input.CalcAverageWidth()){
                
                    std::cout << std::endl << "Calculating average resonance widths ... " << std::endl;
    
                    std::vector<std::pair<double,double>> sWave(xs->excitationEnergies_.size()); 
                    std::vector<std::pair<double,double>> pWave(xs->excitationEnergies_.size());
                    std::vector<std::pair<double,double>> dWave(xs->excitationEnergies_.size());
    
                    
                    ProgressBar pg;
                    pg.start(xs->excitationEnergies_.size());
                    
                    for(unsigned int i = 0; i < xs->excitationEnergies_.size(); i++){
                        pg.update(i);
                        std::cout.precision(3);
                        sWave.at(i) = xs->CalcAverageSWaveResWidth(xs->excitationEnergies_.at(i), 0);
                        pWave.at(i) = xs->CalcAveragePWaveResWidth(xs->excitationEnergies_.at(i), 0);
	                    dWave.at(i) = xs->CalcAverageDWaveResWidth(xs->excitationEnergies_.at(i), 0);
                    }
    
                    pg.update(xs->excitationEnergies_.size());
                    
                    std::cout << std::endl;
                    std::cout << std::endl;
                    std::cout << "Energy\t" << "s width\t" << "s spacing\t" << "p width\t" << "p spacing\t" << "d width\t" << "d spacing\t" << std::endl;
                    std::cout << " [MeV]\t" << "  [meV]\t" << "    [keV]\t" << "  [meV]\t" << "    [keV]\t" << "  [meV]\t" << "    [keV]\t" << std::endl;
    
                    for(unsigned int i = 0; i < xs->excitationEnergies_.size(); i++){
                    
                                std::cout.precision(2);
    
                                std::cout << std::scientific << xs->excitationEnergies_.at(i)
                                                        << "\t" << sWave.at(i).first << "\t" << sWave.at(i).second 
                                                        << "\t" << pWave.at(i).first << "\t" << pWave.at(i).second 
                                                        << "\t" << dWave.at(i).first << "\t" << dWave.at(i).second 
                                                        <<  std::endl;
                            
                    }
                }

                std::cout << std::endl << std::endl << "Calculating for Z: " << it->Z_ << " A: " << it->A_ << " and Particle Type: " << it->pType_ << std::endl;
                
                if(input.CalcXS()) xs->Calculate();
                if(input.CalcXS() && input.PrintXs()) xs->PrintCrossSections();
                if(input.PrintTrans()) xs->PrintTransmissionTerms();
                if(input.CalcRates()) xs->CalculateReactionRates(false);
                if(input.PrintRate() && input.CalcRates()) xs->PrintReactionRates(false);
                if(input.CalcMACS()) xs->CalculateReactionRates(true);
                if(input.PrintMACS() && input.CalcMACS()) xs->PrintReactionRates(true);           
            }
            else
            {
                std::cout << "Could not calculate cross section." << std::endl;    
            }            
            delete xs;
        }
        std::cout << std::fixed;
    }
    

    void RunSingleReaction(const SapphireInput & input){
        std::cout<< std::endl << "Starting calculations for reaction ... " << input.Reaction() << std::endl << std::endl;
        input.SetInputCrossSection();
        input.SetInputDecayer();
        input.SetInputTransitionRate();
        input.SetInputParticleTransmission();
        input.SetInputGammaTransmission();
        input.SetInputLevelDensity();
        GammaTransmissionFunc::SetGnorm(1.0);

        std::vector<int> exitStates(4,-1);
        exitStates[0]=input.g_ExitStates();
        exitStates[1]=input.n_ExitStates();
        exitStates[2]=input.p_ExitStates();
        exitStates[3]=input.a_ExitStates();        

        std::string reactionString(input.Reaction());
        int A = massNumberIntFromString(reactionString);                    
        int Z = atomicNumberIntFromString(reactionString);

        if(!Z && !A){
            std::cout<< std::endl << "Could not get a valid target nucleus from reaction ... "  << std::endl;
            exit(1);
        }

        int pType = pTypeIntFromString(reactionString);
        std::string energyFile = input.EnergyFile();
        bool forRates = input.CalcRates();
        int entranceState = input.EntranceState();
                  

        CrossSection* xs = new CrossSection(Z,A,pType,energyFile,forRates,entranceState,exitStates);
        if(xs->IsValid())
        {
            //xs->CalcEntranceWidth();

            if(input.CalcAverageWidth()){
                
                std::cout << std::endl << "Calculating average resonance widths ... " << std::endl;

                std::vector<std::pair<double,double>> sWave(xs->excitationEnergies_.size()); 
                std::vector<std::pair<double,double>> pWave(xs->excitationEnergies_.size());
                std::vector<std::pair<double,double>> dWave(xs->excitationEnergies_.size());

                
                ProgressBar pg;
                pg.start(xs->excitationEnergies_.size());
                
                for(unsigned int i = 0; i < xs->excitationEnergies_.size(); i++){
                    pg.update(i);
                    std::cout.precision(3);
                    sWave.at(i) = xs->CalcAverageSWaveResWidth(xs->excitationEnergies_.at(i), 0);
                    pWave.at(i) = xs->CalcAveragePWaveResWidth(xs->excitationEnergies_.at(i), 0);
	                dWave.at(i) = xs->CalcAverageDWaveResWidth(xs->excitationEnergies_.at(i), 0);
                }

                pg.update(xs->excitationEnergies_.size());
                
                std::cout << std::endl;
                std::cout << std::endl;
                std::cout << "Energy\t" << "s width\t" << "s spacing\t" << "p width\t" << "p spacing\t" << "d width\t" << "d spacing\t" << std::endl;
                std::cout << " [MeV]\t" << "  [meV]\t" << "    [keV]\t" << "  [meV]\t" << "    [keV]\t" << "  [meV]\t" << "    [keV]\t" << std::endl;

                for(unsigned int i = 0; i < xs->excitationEnergies_.size(); i++){

                            std::cout.precision(2);

                            std::cout << std::scientific << xs->excitationEnergies_.at(i)
                                                    << "\t" << sWave.at(i).first << "\t" << sWave.at(i).second 
                                                    << "\t" << pWave.at(i).first << "\t" << pWave.at(i).second 
                                                    << "\t" << dWave.at(i).first << "\t" << dWave.at(i).second 
                                                    <<  std::endl;
                        
                }

                std::cout << std::fixed;
            }

            if(input.CalcXS()) xs->Calculate();
            if(input.CalcXS() && input.PrintXs()) xs->PrintCrossSections();
            if(input.PrintTrans()) xs->PrintTransmissionTerms();
            if(input.CalcRates()) xs->CalculateReactionRates(false);
            if(input.PrintRate() && input.CalcRates()) xs->PrintReactionRates(false);
            if(input.CalcMACS()) xs->CalculateReactionRates(true);
            if(input.PrintMACS() && input.CalcMACS()) xs->PrintReactionRates(true);  
        }
        else
        {
            std::cout << "Could not calculate cross section." << std::endl;    
        }

        delete xs;
    }

    void RunSingleReaction(std::string &reactionString){
        /** This method will invoke a cross section calculation for a reaction given as reactionString .*/
        
        std::string reaction=reactionString;
        int A = massNumberIntFromString(reactionString);
        int Z = atomicNumberIntFromString(reactionString);
        int pType = pTypeIntFromString(reactionString);
        std::string energyFile;
        bool forRates = false;
        int entranceState = 0;
        std::vector<int> exitStates(4,-1);

        CrossSection* xs = new CrossSection(Z,A,pType,energyFile,forRates,entranceState,exitStates);

        if(xs->IsValid())
        {   
            std::cout << "Calculate cross section for reaction: " << reaction << std::endl;
            xs->Calculate();
            xs->PrintCrossSections();
        }
        else
        {
            std::cout << "Could not calculate cross section." << std::endl;    
        }            
        delete xs;
    }

    void Go(int argc,char *argv[]){
         
        /** At the beginning of the Go() method a clock is started for the measurement of the total calculation time.*/
        auto start = std::chrono::steady_clock::now();
        std::cout << "Module: reaction" << std::endl;
        std::cout << std::endl;
        /** It'll be checked, whether enough cmd line parameters are given, otherwise the printHelp() method is called. */
        if(argc < 3){
            printHelp();
            exit(0);
        }
        /** If enough cmd line parameters are given, then it'll be checked if the inputFile given is an actual file.*/
        if(fexists(argv[2])){              
            /** If the inputFile can be found, a SapphireInput object is created, the inputFile will be printed out
            * via SapphireInput::printInputFile() and the content will be read into the SapphireInput object via
            * SapphireInput::ReadInputFile().
            */
            SapphireInput Input;    
            std::string str(argv[2]);
            Input.PrintIntputFile(str);
            Input.ReadInputFile(str);
            
            /** The next step is to test if a single reaction calculation, or a calculations for reactions given in a
            * reactionFile need to be performed. If the reactionFile is a actual file, the calculation continues with Run().
            * If the reactionFile is not a regular file, then the calculations will be performed with RunSingleReaction().
            * This means, if both keywords are set, the reactionFile option will overwrite the reaction keyword. */

            if(fexists(Input.ReactionFile().c_str()))
            {
                std::cout << std::endl << "Starting calculations for reactions in file ... " << Input.ReactionFile().c_str() << std::endl<<std::endl;
                Input.Reaction("");
                Input.PrintIntputParameters("CrossSection");
                Run(Input);
            }
            else
            {
                std::cout << std::endl << "No valid reactionFile was given...\nUsing " << Input.Reaction() <<"..." << std::endl;
                Input.PrintIntputParameters("CrossSection");
                RunSingleReaction(Input);
            }          
                     
        }
        else{
            /** If no regular inputFile is given in the first place, the code assumes that the cmd line parameter is a
            * reactionString. In this case RunSingleReaction will be called by passing the cmd line parameter.*/
            std::string reactionString(argv[2]);
            RunSingleReaction(reactionString);
        }

        /** At the end of the Go() method, the clock is stopped and the total calculation time is printed to std::cout.
        * 
        * ---
        */
        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << std::endl << "Total calculation time: " << elapsed_seconds.count() << "s\n";
    } 
}
