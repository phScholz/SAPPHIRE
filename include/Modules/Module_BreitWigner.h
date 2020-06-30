/**
 * @file Module_BreitWigner.h
 * @brief Contains methods to calculate HF and BW cross sections from the input of beta-delayed neutron experiments
 */

#pragma once
#include <iostream>
#include <vector>
#include <cmath>

#include "WidthEntry.h"
#include "CompoundState.h"

namespace Module_BreitWigner{
    /** 
    *   @brief Top level function to call from main for the reaction Module.
    *   @param argc Number of cmd line parameters of Sapphire.
    *   @param argv The array which contains the cmd line parameters of Sapphire.
    */
    void Go(int argc,char *argv[]); 

    /**
     * @brief Function to calculate the compound decay widths from a partial width file
     * 
     */
    void CalcWidths(std::string reaction, std::string file);

    /**
     * @brief Calculate BW cross section if energyFile is given
     * 
     * @param reaction 
     * @param widths 
     * @param energies 
     * 
     * @todo 
     * - Separate this method from printing to std::out
     * - Create a write method for storing results in a file
     * - formatting of output
     * - create a class as container for the results
     * - simultanously calculate HF cross section
     * - adjusting level density
     * - using input partial width in total width
     */
    void CalcBreitWigner(std::string reaction, std::string widths, std::string energies);

    double BW(double G1, double G2, double Gtot, double density, double E1, double E);

    /**
    * @brief Check wheter a string represents an actual file.
    * @param filename String with the supposedly path to a file.
    * @return True if the file exists; False if it doesn't.
    */
    bool fexists(const char *filename);
    
    /**
    * @brief Print the information of Compound States to std::cout.
    * 
    * @param compound A CompoundStates object.
    */
    void PrintCompoundStates(std::vector<CompoundState> states);

    /**
     * @brief Read the partial width file and return a vector of WidthEntry objects for each line.
     * 
     * @param file 
     * @return std::vector<WidthEntry> Vector of WidthEntry objects
     */
    std::vector<WidthEntry> ReadWidthFile(std::string file);

    /**
     * @brief Read the energy file and return a vector of energy doubles for each line.
     * 
     * @param file 
     * @return std::vector<double> of energies
     */
    std::vector<double> ReadEnergyFile(std::string file);

    /**
     * @brief Printing help information to std::cout
     */
    void PrintHelp();

    /**
     * @brief Get a projectile integer from a projectile string
     * @param reactionString The part of the reactionString which defines the projectile.
     * @return Integer which represents the pType: 0 = gamma, 1 = neutron, 2 = proton, 3 = alpha.
     */
    int pTypeIntFromString(std::string &reactionString);

        /**
     * @brief Get the pTypeString from the reactionString.
     * @param reactionString A string from the cmd line which represents the reaction, e.g., 60Fe+p
     * @returns The pType part of the reactionString, e.g. "p"
     */
    std::string pTypeStringFromString(std::string &reactionString);

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


}