/**
 * @file Module_NLD.h
 * @brief Contains the declarations for the NLD module
 */

#pragma once
#include <string>
#include "LevelDensity/LevelDensity.h"

namespace Module_NLD{
    /** 
    *   @brief Top level function to call from main for the reaction Module.
    *   @param argc Number of cmd line parameters of Sapphire.
    *   @param argv The array which contains the cmd line parameters of Sapphire.
    */
    void Go(int argc,char *argv[]); 

    /**
     * @brief Function to get NLD
     */
    void GetNLD(std::string isotope, int nldmodel);

    /**
     * @brief Print NLD for model 0
     */
    template<class T>
    void Print(std::string isotope);

    /**
     * @brief Printing help information to std::cout
     */
    void PrintHelp();

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