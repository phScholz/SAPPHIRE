/**
 * @file Module_CrossSection.h
 * @brief Reimplementation of the cross section calculations routine of Sapphire as a module.
 * @date 2020
 * @author Philipp Scholz, <pscholz@outlook.de>
 * 
 */
#pragma once
#include <vector>
#include <map>
#include <string>
#include "SapphireInput.h"

namespace Module_CrossSection{
    typedef struct EntrancePairs {
        EntrancePairs(int Z,int A,int pType) {
            Z_=Z;
            A_=A;
            pType_=pType;
        };
        int Z_;
        int A_;
        int pType_;
    } EntrancePairs;

    /** 
    *   @brief Top level function to call from main for the reaction Module.
    *   @param argc Number of cmd line parameters of Sapphire.
    *   @param argv The array which contains the cmd line parameters of Sapphire.
    */
    void Go(int argc,char *argv[]); 

    /**
    *   @brief Function to perform cross section calculations on the basis of a SapphireInput object.
    *   @param input This is a SapphireInput object which should contain all parameters for the calculation of the cross section.
    */
    void Run(const SapphireInput & input); 

    /**
    *   @brief Function to perform a single reaction cross section calculations on the basis of a SapphireInput object.
    *   @param input This is a SapphireInput object which should contain all parameters for the calculation of the cross section.
    */
    void RunSingleReaction(const SapphireInput & input);

    /**
    *   @brief Function to perform a single reaction cross section calculations on the basis of a reactionString with default parameters.
    *   @param reactionString This is a std::string object which contains the reaction in a format similar to "60Fe+p".
    */
    void RunSingleReaction(std::string reactionString);

    /**
     * @brief Method to read in the Entrance Pairs given in the reactionFile
     * @param entrancePairs A std::vector object which contains entrancePairs
     * @param reactionFile The file which contains a list of target charge and mass, as well as projectile type.
     */
    void readEntrancePairs(std::vector<EntrancePairs> & entrancePairs, std::string reactionFile);

    /**
    * @brief method to print the Entrance Pairsto stdout
    */
    void PrintEntrancePairs(std::vector<EntrancePairs> & entrancePairs);

    /**
     * @brief Check wheter a string represents an actual file.
     * @param filename String with the supposedly path to a file.
     * @return True if the file exists; False if it doesn't.
     */
    bool fexists(const char *filename);
    
    /**
     * @brief Get the massNumberString from the reactionString.
     * @param reactionString A string from the cmd line which represents the reaction, e.g., 60Fe+p
     * @returns The mass number part of the reactionString, e.g. "60"
     */
    std::string massNumberStringFromString(std::string reactionString);

    /**
     * @brief Get the massNumberInt from the reactionString.
     * @param reactionString A string from the cmd line which represents the reaction, e.g., 60Fe+p
     * @returns The mass number part of the reactionString, e.g. 60, if successful. Otherwise returns 0.
     */
    int massNumberIntFromString(std::string reactionString);

    /**
     * @brief Get the pTypeString from the reactionString.
     * @param reactionString A string from the cmd line which represents the reaction, e.g., 60Fe+p
     * @returns The pType part of the reactionString, e.g. "p"
     */
    std::string pTypeStringFromString(std::string reactionString);

    /**
     * @brief Get a projectile integer from a projectile string
     * @param reactionString The part of the reactionString which defines the projectile.
     * @return Integer which represents the pType: 0 = gamma, 1 = neutron, 2 = proton, 3 = alpha.
     */
    int pTypeIntFromString(std::string reactionString);

    /**
     * @brief Get the atomic number as string from reactionString.
     * @param reactionString A string from the cmd line which represents the reaction, e.g., 60Fe+p
     * @returns The atomic number part of the reactionString, e.g. Fe
     */
    std::string atomicNumberStringFromString(std::string reactionString);

    /**
     * @brief Get the atomic number as int from reactionString.
     * @param reactionString A string from the cmd line which represents the reaction, e.g., 60Fe+p
     * @returns The atomic number part of the reactionString, e.g. 28 other wise returns 0.
     */
    int atomicNumberIntFromString(std::string reactionString);

    /**
    *   @brief Funtion to print out the help information for the reaction module.
    */
    void printHelp();
}