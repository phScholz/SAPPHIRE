#include "Module_CrossSection.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>

/*Time measurement*/
#include <chrono>
#include <ctime>

#include "NuclearMass.h"
#include "CrossSection.h"
#include "SapphireInput.h"
#include "Decayer.h"
#include "TransitionRateFunc.h"
#include "ParticleTransmissionFunc.h"
#include "GammaTransmissionFunc.h"



namespace Module_CrossSection{
 

    bool fexists(const char *filename) {
        std::ifstream ifile(filename);
    return (bool)ifile;
    }

    std::string pTypeStringFromString(std::string reactionString){
        std::string massNumberString=massNumberStringFromString(reactionString);
        std::string atomicNumberString=atomicNumberStringFromString(reactionString);

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
        std::string atomicNumberString = atomicNumberStringFromString(reactionString);
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
        std::cout << "\tAX+y           - reaction string, e.g. 60Fe+p (not yet implemented)" << std::endl;
        std::cout << "\tInputFile      - InputFile (not yet implemented)" << std::endl;
        std::cout << std::endl;
    }

    void readEntrancePairs(std::vector<EntrancePairs> & entrancePairs, std::string reactionFile){
        std::ifstream in(reactionFile.c_str());
        int Z,A,pType =0;
        if(!in) {
	        std::cout << "Could not open " << reactionFile << " for reading." << std::endl;
	        exit(1);
        } else {
	        std::cout << "Reading nuclei from " << reactionFile << "." << std::endl;

	        while(!in.eof()) {
	            std::string line;
	            std::getline(in,line);
	            
                if(!in.eof()) {
	                std::istringstream stm(line);
	                if(stm >> Z >> A >> pType)
	                    entrancePairs.push_back(EntrancePairs(Z,A,pType));
	            }
	        }
	                
        in.close();
        }
    }

    void Run(SapphireInput * input){
        //Copy the input parameters to respective classes
        //I know, this is still not really intuitive ... but I am working on it.
        CrossSection::SetResidualAlpha(input->ResidualAlpha());
        CrossSection::SetResidualProton(input->ResidualProton());
        CrossSection::SetResidualNeutron(input->ResidualNeutron());
        CrossSection::SetResidualGamma(input->ResidualGamma());
        CrossSection::SetCalculateGammaCutoff(input->CalculateGammaCutoff());
        TransitionRateFunc::SetGammaCutoffEnergy(input->g_CutoffEnergy());
        Decayer::SetCrossSection(true);
        Decayer::SetMaxL(input->DecayerMaxL());
        
        std::vector<int> exitStates(4,-1);
        exitStates[0]=input->g_ExitStates();
        exitStates[1]=input->n_ExitStates();
        exitStates[2]=input->p_ExitStates();
        exitStates[3]=input->a_ExitStates();
        
        ParticleTransmissionFunc::SetAlphaFormalism(input->a_Formalism());
        ParticleTransmissionFunc::SetProtonFormalism(input->p_Formalism());
        ParticleTransmissionFunc::SetNeutronFormalism(input->n_Formalism());
        ParticleTransmissionFunc::SetPorterThomas(input->PorterThomas_p());

        GammaTransmissionFunc::SetEGDRType(input->g_Formalism());
        GammaTransmissionFunc::SetPorterThomas(input->PorterThomas_g());

        std::vector<EntrancePairs> entrancePairs;
        readEntrancePairs(entrancePairs,input->ReactionFile());

        for(auto it = std::begin(entrancePairs); it != std::end(entrancePairs); ++it){}
            CrossSection* xs = new CrossSection(it.Z_,it.A_, it.pType_,input->EnergyFile(),input->CalcRates());
            if(xs->IsValid())
            {
                std::cout << "Calculating for Z: " << it[i].Z_ << " A: " << it[i].A_ << " and Particle Type: " << it[i].pType_ << std::endl;
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

    void RunSingleReaction(SapphireInput * input){
        std::cout<< std::endl << "Starting calculations for reaction ... " << input->Reaction() << std::endl;
            int A = massNumberIntFromString(input->Reaction());
            int Z = atomicNumberIntFromString(input->Reaction());

            if(!Z && !A){
                std::cout<< std::endl << "Could not get a valid target nucleus from reaction ... "  << std::endl;
                exit(1);
            }

            int pType = pTypeIntFromString(input->Reaction());
            std::string energyFile = input->EnergyFile();
            bool forRates = input->CalcRates();
            int entranceState = input->EntranceState();
            std::vector<int> exitStates(4,-1);
            exitStates[0]=input->g_ExitStates();
            exitStates[1]=input->n_ExitStates();
            exitStates[2]=input->p_ExitStates();
            exitStates[3]=input->a_ExitStates();

            std::cout << "Input Values For Cross Section:"   << std::endl
		            << std::setw(14) << "Z:"               << std::setw(12) 
		            << Z      << std::setw(0) << std::endl
		            << std::setw(14) << "A:"               << std::setw(12) 
		            << A      << std::setw(0) << std::endl;
                
                
            if(pType==0) 
	            std::cout << std::setw(14) << "projectile:"  << std::setw(12) << "g" << std::setw(0) << std::endl;
            else if(pType==1) 
	            std::cout << std::setw(14) << "projectile:"  << std::setw(12) << "n" << std::setw(0) << std::endl;
            else if(pType==2) 
	            std::cout << std::setw(14) << "projectile:"  << std::setw(12) << "p" << std::setw(0) << std::endl;
            else if(pType==3) 
	            std::cout << std::setw(14) << "projectile:"  << std::setw(12) << "a" << std::setw(0) << std::endl;
                
            std::cout << "Starting Cross Section Calculation..." << std::endl;           

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

    void RunSingleReaction(std::string reactionString){
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
        auto start = std::chrono::steady_clock::now();
        std::cout << std::endl;
        std::cout << "Sapphire reaction" << std::endl;
        std::cout << "*****************" << std::endl;

        if(argc < 3){
            printHelp();
            exit(0);
        }
        
        if(fexists(argv[2])){              
            SapphireInput* Input = new SapphireInput();            
            Input->printIntputFile();
            Input->ReadInputFile(str(argv[2]));
            Input->printIntputParameters();

            /*Defined in Setup.cpp ... Should not be in another file*/
            //ElementTable NuclearMass::elementTable_; 
            //MassTable NuclearMass::massTable_;
            //GDRTable GammaTransmissionFunc::gdrTable_;
            //LevelsTable NuclearLevels::levelsTable_;

            if(fexists(Input->ReactionFile().c_str()))
            {
                std::cout << std::endl << "Starting calculations for reactions in file ... " << argv[2] << std::endl;
                Run(Input);
            }
            else
            {
                std::cout << std::endl << "No valid reactionFile was given...\nUsing " << Input->Reaction <<"..." << std::endl;
                RunSingleReaction(Input);
            }          
            delete Input;           
        }
        else{
            RunSingleReaction(argv[2]);
        }

        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << std::endl << "elapsed time: " << elapsed_seconds.count() << "s\n";
    } 
}
