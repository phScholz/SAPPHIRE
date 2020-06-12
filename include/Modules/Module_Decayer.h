/**
 * @file Module_Decayer/Decayer.h
 * @brief Reimplementation of the decay routine of Sapphire as a module.
 * @date 2020
 * @author Philipp Scholz, <pscholz@outlook.de>
 * 
 */

#pragma once
#include <vector>
#include <map>
#include <string>
#include "SapphireMPITypes.h"
#include "SapphireInput.h"
#include "SPDistribution.h"

namespace Module_Decayer{
    void Go(int argc,char *argv[]); /**< Top level function to call from main*/
    void RunSingle(const SapphireInput & input); /**< Declaration of the main function of the Decayer Module*/
    void RunDist(const SapphireInput & input); /**< Declaration of the main function with spin distribution of the Decayer Module*/

    /**
     * @brief Reading Spin-Parity Distribution from a File
     * @param file Path to the spin-parity distribution file
     * @param dist SPDistribution paramter by reference.
     * @return True by success; False by failor.
     */
    bool ReadDist(std::string file, SPDistribution &dist);

    /**
     * @brief Check wheter a string represents an actual file.
     * @param filename String whith the supposedly path to a file.
     * @return True if the file exists; False if it doesn't.
     */
    bool fexists(const char *filename);

    /**
     * @brief Get the atomic number as string from isotopeString.
     * @param isotopeString A string from the cmd line which represents the reaction, e.g., 60Fe+p
     * @returns The atomic number part of the isotopeString, e.g. Fe
     */
    std::string atomicNumberStringFromString(std::string &isotopeString);

    /**
     * @brief Get the atomic number as int from isotopeString.
     * @param isotopeString A string from the cmd line which represents the reaction, e.g., 60Fe+p
     * @returns The atomic number part of the isotopeString, e.g. 28 other wise returns 0.
     */
    int atomicNumberIntFromString(std::string &isotopeString);

    /**
     * @brief Get the massNumberString from the isotopeString.
     * @param isotopeString A string from the cmd line which represents the reaction, e.g., 60Fe+p
     * @returns The mass number part of the isotopeString, e.g. "60"
     */
    std::string massNumberStringFromString(std::string &isotopeString);

    /**
     * @brief Get the massNumberInt from the isotopeString.
     * @param isotopeString A string from the cmd line which represents the reaction, e.g., 60Fe+p
     * @returns The mass number part of the isotopeString, e.g. 60, if successful. Otherwise returns 0.
     */
    int massNumberIntFromString(std::string &isotopeString);

    //void masterProcess(boost::mpi::communicator& world,InitialNucleusData initalNucleus, int suffixNo,int events);
    //void slaveProcess(boost::mpi::communicator& world,InitialNucleusData initialNucleus);

    void printHelp();
}