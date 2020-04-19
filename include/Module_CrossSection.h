/**
 * @file Module_CrossSection.h
 * @brief Reimplementation of the cross section calculations routine of Sapphire as a module.
 * @date 2020
 * @author Philipp Scholz, <pscholz@outlook.de>
 * 
 */

#ifndef MODULE_CROSSSECTION_H
#define MODULE_CROSSSECTION_H
#endif
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

    void Go(int argc,char *argv[]); /**< Top level function to call from main*/
    void Run(const SapphireInput & input); /**< Declaration of the main function of the CrossSection Module*/

    void RunSingleReaction(const SapphireInput & input);
    void RunSingleReaction(std::string reactionfile);

    /**
     * @brief method to read in the Entrance Pairs given in the reactionFile
     */
    void readEntrancePairs(std::vector<EntrancePairs> & entrancePairs, std::string reactionFile);

    /**
    * @brief method to print the Entrance Pairsto stdout
    */
    void PrintEntrancePairs(std::vector<EntrancePairs> & entrancePairs);

    /**
     * @brief Check wheter a string represents an actual file.
     * @param filename String whith the supposedly path to a file.
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
     * @param projectileString The part of the reactionString which defines the projectile.
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

    void printHelp();
}