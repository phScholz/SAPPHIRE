/**
 * @file Module_Compound.h
 * @brief Contains the declarations for the Compound module to calculate decay widths for compound states.
 */

#pragma once
#include <string>
#include "LevelDensity/LevelDensity.h"
#include "CompoundStates.h"

namespace Module_Compound{
    /** 
    *   @brief Top level function to call from main for the reaction Module.
    *   @param argc Number of cmd line parameters of Sapphire.
    *   @param argv The array which contains the cmd line parameters of Sapphire.
    */
    void Go(int argc,char *argv[]); 

    /**
     * @brief Print the information of Compound States to std::cout.
     * 
     * @param compound A CompoundStates object.
     */
    void PrintCompoundStates(std::vector<CompoundState> states);

    /**
     * @brief Write the information of Compound States to a file.
     * 
     * @param compound A CompoundStates object.
     * @param file The path to the output file
     */
    void WriteCompoundStates(CompoundStates * compound, std::string file);

    /**
    * @brief Check wheter a string represents an actual file.
    * @param filename String with the supposedly path to a file.
    * @return True if the file exists; False if it doesn't.
    */
    bool fexists(const char *filename);

    /**
     * @brief Printing help information to std::cout
     */
    void PrintHelp();

    /**
     * @brief Function to calculate the compound decay widths
     * 
     */
    void CalcWidths(std::string reaction, double energy);


    /**
     * @brief Same as CalcWidths() but for only one compound state
     * 
     * @param reactionString Capture reaction
     * @param energy energy file or single energy
     * @param state spin and parity e.g. "1/2+"
     */
    void CalcWidths(std::string reactionString, double energy, std::string state);

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

    double spinDoubleFromString(std::string &jpi);

    int parityIntFromString(std::string &jpi);

    void parseJPiString(double & spin, int & parity, std::string & jPiString);
}