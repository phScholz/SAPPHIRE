/**
 * @file Module_OldSapphire.h
 * @brief Reimplementation of the old Sapphire main routine as a module in the new version.
 * @date 2020
 * @author Philipp Scholz, <pscholz@outlook.de>
 * 
 */

#pragma once
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <algorithm>
#include "DecayController.h"
#include "NuclearMass.h"
#include "DecayResults.h"
#include "CrossSection.h"
#include "omp.h" /** Currently only used for the Decayer*/
#include "SapphireMPITypes.h"
#include "TransitionRateFunc.h"
#include "ParticleTransmissionFunc.h"
#include "GammaStrength/GammaTransmissionFunc.h"

namespace Module_OldSapphire{
    void Go(int argc,char *argv[]); /**< Top level function to call from main*/
    int oldSapphire_MPI(int argc,char *argv[]); /**< Top level function to call from main*/
    int oldSapphire(int argc,char *argv[]); /**< Top level function to call from main*/

    void parseCommandLineForOptions(std::vector<std::string>& args,
                                        int& suffixNo,
                                        bool &preEq,
				                        int& numPiParticles,
                                        int& numPiHoles,
                                        int& numNuParticles,
                                        int& numNuHoles);


    /**
    * @brief CMD line parameters are parsed for the Cross section module
    * @param args cmd line string
    * @param Z nuclear charge number
    * @param A nuclear mass number
    * @param pType Projectile
    * @param energyFile Reference to the string for the path to the EnergyFile
    * @param asciiIn Boolean which shows if there is a asciiFile or not
    * @returns True or False
    */
    bool parseCommandLineForXS(std::vector<std::string>& args,int& Z, int&A, int& pType, std::string& energyFile, bool asciiIn);

    
    /**
    * @brief CMD line parameters are parsed for the Decay module
    * @param args cmd line string
    * @param Z nuclear charge number
    * @param A nuclear mass number
    * @param J Reference to a spin double
    * @param Pi Reference to a parity int
    * @param lowEnergy Reference to a lowEnergy double
    * @param highEnergy  Reference to a highEnergy double
    * @param events  Reference to the events int.
    * @returns True or False
    */
    bool parseCommandLineForDecay(std::vector<std::string>& args, int& Z, int& A, double& J, int& Pi, double& lowEnergy, double& highEnergy, int& events);

    
//    void masterProcess(boost::mpi::communicator& world,InitialNucleusData initalNucleus, int suffixNo,int events);
  
//    void slaveProcess(boost::mpi::communicator& world,InitialNucleusData initialNucleus);
    
    /**
     * @brief Print help for the OldSapphire Module to std::cout.
     */
    void printHelp(); /**< The old printHelp() function*/
}