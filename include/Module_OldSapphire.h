/**
 * @file Module_OldSapphire.h
 * @brief Reimplementation of the old Sapphire main routine as a module in the new version.
 * @date 2020
 * @author Philipp Scholz, <pscholz@outlook.de>
 * 
 */

#ifndef MODULE_OLDSAPPHIRE_H
#define MODULE_OLDSAPPHIRE_H
#endif
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
#ifndef MPI_BUILD
#include "CrossSection.h"
#include "omp.h" /** Currently only used for the Decayer*/
#else
#include "SapphireMPITypes.h"
#include <boost/mpi.hpp>
#endif
#include "TransitionRateFunc.h"
#include "ParticleTransmissionFunc.h"
#include "GammaTransmissionFunc.h"

namespace Module_OldSapphire{
    void Go(int argc,char *argv[]); /**< Top level function to call from main*/
    int oldSapphire_MPI(); /**< Top level function to call from main*/
    int oldSapphire(); /**< Top level function to call from main*/

    #ifdef MPI_BUILD
    void parseCommandLineForOptions(std::vector<std::string>& args,
                                        int& suffixNo,
                                        bool &preEq,
				                        int& numPiParticles,
                                        int& numPiHoles,
                                        int& numNuParticles,
                                        int& numNuHoles);
    #endif

    #ifndef MPI_BUILD
    void parseCommandLineForOptions(std::vector<std::string>& args,
                                    int& suffixNo,
                                    bool &preEq,
				                    int& numPiParticles,
                                    int& numPiHoles,
                                    int& numNuParticles,
                                    int& numNuHoles,
                                    bool& calcAverageWidth,
                                    bool& calcRates, 
                                    bool& asciiIn,
				                    std::string& inFile,
                                    int& entranceState,
                                    std::vector<int>& exitStates,
				                    bool& printTrans);
    #endif

    #ifndef MPI_BUILD
    /**
    * @brief CMD line parameters are parsed for the Cross section module
    * @param args cmd line string
    * @param Z nuclear charge number
    * @param A nuclear mass number
    * @param pType Projectile
    * @param energyFile Reference to the string for the path to the EnergyFile
    * @param asciiIn Boolean which shows if there is a asciiFile or not
    * @param highEnergy  Reference to a highEnergy double
    * @param events  Reference to the events int.
    * 
    * @returns True or False
    * 
    * @note This has to move to another class in future. This does not belong in a main file but in the class file of the respective module.
    * 
    */
    bool parseCommandLineForXS(std::vector<std::string>& args,int& Z, int&A, 
			   int& pType, std::string& energyFile, bool asciiIn);
    #endif

    
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
     * 
     * @returns True or False
     * 
     * @note This has to move to another class in future. This does not belong in a main file but in the class file of the respective module.
     * 
     */
    bool parseCommandLineForDecay(std::vector<std::string>& args, 
			      int& Z, int& A, double& J, int& Pi, 
			      double& lowEnergy, double& highEnergy,
			      int& events)

    #ifdef MPI_BUILD
    void masterProcess(boost::mpi::communicator& world,InitialNucleusData initalNucleus, int suffixNo,int events);
    #endif

    #ifdef MPI_BUILD
    void slaveProcess(boost::mpi::communicator& world,InitialNucleusData initialNucleus);
    #endif

    void printHelp(); /**< The old printHelp() function*/
}