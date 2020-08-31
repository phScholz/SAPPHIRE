#include "Modules/Module_BreitWigner.h"
#include "Databases/NuclearMass.h"
#include "CrossSection.h"
#include "CompoundStates.h"
#include "CompoundState.h"
#include "WidthEntry.h"
#include "Constants.h"

#include <iostream>
#include <iomanip> 
#include <string>
#include <fstream>
#include <sstream>
#include <chrono>
#include <cmath>
#include <vector>

namespace Module_BreitWigner{
    void Go(int argc, char *argv[]){
        auto start = std::chrono::steady_clock::now();
        std::cout << std::endl << "Module: breitWigner" << std::endl;
        std::cout << std::endl;
        if(argc < 3){
            std::cout << std::endl << "Not enough input parameters ..." << std::endl;
            PrintHelp();
            exit(0);
        }

        if(argc == 4){
            std::string reaction(argv[2]);
            std::string file(argv[3]);
            CalcWidths(reaction, file);
        }

        if(argc == 5){
            std::string reaction(argv[2]);
            std::string widthFile(argv[3]);
            std::string energyFile(argv[4]);
            CalcBreitWigner(reaction, widthFile, energyFile);
        }

        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << std::endl << "Total calculation time: " << elapsed_seconds.count() << "s\n";


    }

    void CalcBreitWigner(std::string reactionString, std::string widthFile, std::string energyFile){
        int A = massNumberIntFromString(reactionString);
        int Z = atomicNumberIntFromString(reactionString);
        int pType = pTypeIntFromString(reactionString);

        std::vector<double> energies = ReadEnergyFile(energyFile);
        std::vector<WidthEntry> widths = ReadWidthFile(widthFile);

        for(auto it = widths.begin(); it != widths.end(); ++it){
            CompoundStates compound;
            std::ofstream fp;
            fp.open("output/energyFile.dat");
            fp << it->energy_ << std::endl;
            fp.close();
            
            energyFile="output/energyFile.dat";

            bool forRates = false;
            int entranceState = 0;
            std::vector<int> exitStates(4,-1);

            std::cout << std::endl <<   "Calculating Compound states for ..." << std::endl 
                                   <<   "         A: " << A << std::endl 
                                   <<   "         Z: " << Z << std::endl
                                   <<   "Projectile: " << pType << std::endl
                                   <<   "     State: " << it->jPiString_ << std::endl
                                   <<   "    Energy: " << it->energy_ << std::endl;

            CrossSection* xs = new CrossSection(Z,A,pType,energyFile,forRates,entranceState,exitStates);

            if(xs->IsValid())
            {
                compound = xs->CalcCompoundWidth(it->spin_, it->parity_);
                //Setting the neutron groundstate width to the experimental value.
                double oldGSWidth = compound.states.at(0).GroundStateWidth("neutron");
                double oldTotalWidth = compound.states.at(0).TotalWidth("neutron");
                compound.states.at(0).GroundStateWidth("neutron",it->width_);
                compound.states.at(0).TotalWidth("neutron",oldTotalWidth - oldGSWidth + it->width_);

                PrintCompoundStates(compound.states);
                
                std::cout << std::endl;
                std::cout << "\t"
                          << std::setw(6) << "Energy [MeV]" << "\t"
                          << std::setw(6) << "Crosssection [b]" << "\t"
                          <<std::endl;

                

                for(auto erg = energies.begin(); erg != energies.end(); ++erg){
                    double bw = BW(it->width_, compound.states.at(0).TotalWidth("gamma"), compound.states.at(0).TotalWidth("all"), compound.states.at(0).Density(), it->energy_, *erg);
                    double sigma = xs->PreFactor() / *erg * (2*it->spin_+1) * bw;
                    std::cout << "\t"
                              << std::setw(6) << std::scientific << *erg << "\t"
                              << std::setw(6) << std::scientific << sigma
                              << std::fixed << std::endl;
                }
            }
        }
    }

    double BW(double G1, double G2, double Gtot, double density, double E1, double E){
        return G1*G2 / ( 4*pi*density*density * (E1-E) * (E1-E) + Gtot*Gtot/4.);
    }

    void CalcWidths(std::string reactionString, std::string file){        
        int A = massNumberIntFromString(reactionString);
        int Z = atomicNumberIntFromString(reactionString);
        int pType = pTypeIntFromString(reactionString);
        std::string energyFile;

        //Read Partial width file
        std::vector<WidthEntry> Input = ReadWidthFile(file);
        
        //For each entry in the partial width file, get compound states
        for(auto it = Input.begin(); it != Input.end(); ++it){
            
            std::ofstream fp;
            fp.open("output/energyFile.dat");
            fp << it->energy_ << std::endl;
            fp.close();
            
            energyFile="output/energyFile.dat";

            bool forRates = false;
            int entranceState = 0;
            std::vector<int> exitStates(4,-1);

            std::cout << std::endl <<   "Calculating Compound states for ..." << std::endl 
                                   <<   "         A: " << A << std::endl 
                                   <<   "         Z: " << Z << std::endl
                                   <<   "Projectile: " << pType << std::endl
                                   <<   "     State: " << it->jPiString_ << std::endl
                                   <<   "    Energy: " << it->energy_ << std::endl;

            CrossSection* xs = new CrossSection(Z,A,pType,energyFile,forRates,entranceState,exitStates);

            if(xs->IsValid())
            {
                CompoundStates compound = xs->CalcCompoundWidth(it->spin_, it->parity_);
                PrintCompoundStates(compound.states);
            }

        }
    }

    bool fexists(const char *filename){
        std::ifstream ifile(filename);
        return (bool)ifile;
    }

    void PrintCompoundStates(std::vector<CompoundState> states){
        std::cout     << std::setw(6) <<  std::fixed << "Mass" << "\t"
                      << std::setw(6) <<  std::fixed << "Charge" << "\t"
                      << std::setw(6) <<  std::fixed << "Spin" << "\t"
                      << std::setw(6) <<  std::fixed << "Parity" << "\t"
                      << std::setw(6) <<  std::fixed << "Energy" << "\t"
                      << std::setw(9) <<  std::scientific << "Density" << "\t"
                      << std::setw(6) <<  std::scientific << "g GS Width  " << "\t"
                      << std::setw(6) <<  std::scientific << "n GS Width  " << "\t"
                      << std::setw(6) <<  std::scientific << "p GS Width  " << "\t"
                      << std::setw(6) <<  std::scientific << "a GS Width  " << "\t"
                      << std::setw(6) <<  std::scientific << "g TotalWidth" << "\t"
                      << std::setw(6) <<  std::scientific << "n TotalWidth" << "\t"
                      << std::setw(6) <<  std::scientific << "p TotalWidth" << "\t"
                      << std::setw(6) <<  std::scientific << "a TotalWidth" << "\t"                                            
                      << std::setw(6) <<  std::scientific << "TotalWidth" << "\t"
                      << std::fixed << std::endl;
        
        std::cout     << std::setw(6) <<  std::fixed << "" << "\t"
                      << std::setw(6) <<  std::fixed << "" << "\t"
                      << std::setw(6) <<  std::fixed << "" << "\t"
                      << std::setw(6) <<  std::fixed << "" << "\t"
                      << std::setw(6) <<  std::fixed << "[MeV]" << "\t"
                      << std::setw(9) <<  std::scientific << "[1/MeV]" << "\t"
                      << std::setw(6) <<  std::scientific << "[MeV]" << "\t\t"
                      << std::setw(6) <<  std::scientific << "[MeV]" << "\t\t"
                      << std::setw(6) <<  std::scientific << "[MeV]" << "\t\t"
                      << std::setw(6) <<  std::scientific << "[MeV]" << "\t\t"
                      << std::setw(6) <<  std::scientific << "[MeV]" << "\t\t"
                      << std::setw(6) <<  std::scientific << "[MeV]" << "\t\t"
                      << std::setw(6) <<  std::scientific << "[MeV]" << "\t\t"
                      << std::setw(6) <<  std::scientific << "[MeV]" << "\t\t"                                            
                      << std::setw(6) <<  std::scientific << "[MeV]" << "\t\t"
                      << std::fixed << std::endl;

        for(auto it = states.begin(); it != states.end(); ++it){
            std::cout << std::setw(6) <<  std::fixed << it->A() << "\t"
                      << std::setw(6) <<  std::fixed << it->Z() << "\t"
                      << std::setw(6) <<  std::fixed << it->J() << "\t"
                      << std::setw(6) <<  std::fixed << ((it->P()>0)? "+" :"-" ) << "\t"
                      << std::setw(6) <<  std::fixed << it->E() << "\t"
                      << std::setw(9) <<  std::scientific << it->Density() << "\t"
                      << std::setw(6) <<  std::scientific << it->GroundStateWidth("gamma") << "\t"
                      << std::setw(6) <<  std::scientific << it->GroundStateWidth("neutron") << "\t"
                      << std::setw(6) <<  std::scientific << it->GroundStateWidth("proton") << "\t"
                      << std::setw(6) <<  std::scientific << it->GroundStateWidth("alpha") << "\t"
                      << std::setw(6) <<  std::scientific << it->TotalWidth("gamma") << "\t"
                      << std::setw(6) <<  std::scientific << it->TotalWidth("neutron") << "\t"
                      << std::setw(6) <<  std::scientific << it->TotalWidth("proton") << "\t"
                      << std::setw(6) <<  std::scientific << it->TotalWidth("alpha") << "\t"                       
                      << std::setw(6) <<  std::scientific << it->TotalWidth("all") << "\t"     
                      << std::fixed << std::endl;
        }
    }

    std::vector<WidthEntry> ReadWidthFile(std::string file){
        std::vector<WidthEntry> content;        

        std::cout << "Reading partial width file ... " << file << std::endl;

        std::ifstream in(file.c_str());
        
        if(!in){
            std::cout << "Could not read partial width file." << std::endl;
            exit(1);
        }

        std::string line;

        while(!in.eof()){
            std::getline(in, line);
            WidthEntry dummy(line);
            content.push_back(dummy);
        }

        return content;
    }

    std::vector<double> ReadEnergyFile(std::string file){
        std::vector<double> energies;

        std::cout << "Reading energy file ... " << file << std::endl;

        std::ifstream in(file.c_str());
        
        if(!in){
            std::cout << "Could not read energy file." << std::endl;
            exit(1);
        }

        std::string line;

        while(!in.eof()){
            std::getline(in, line);
            std::istringstream lineStream(line);
            double energy;
            lineStream >> energy;
            energies.push_back(energy);
        }

        return energies;
    }

    void PrintHelp(){
        std::cout << std::endl << "\tSyntax:        sapphire breitWigner <xxY+z> <partialWidthsFile> <energyFile>" << std::endl;        
        std::cout << std::endl << "\t               sapphire breitWigner <xxY+z> <partialWidthsFile>" << std::endl;
        std::cout << std::endl << "\t               sapphire breitWigner <InputFile>" << std::endl;
	    std::cout << std::endl;
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

}