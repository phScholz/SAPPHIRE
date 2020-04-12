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

/* Includes from Setup.cpp*/
#include "NuclearMass.h"
#include "GammaTransmissionFunc.h"
#include "NuclearLevels.h"
#include "Decayer.h"
#include "Sapphire_config.h"
#include "TransitionRateFunc.h"
#ifndef MPI_BUILD
#include "CrossSection.h"
#endif
#include "PreEqDecayer.h"
#include "ParticleTransmissionFunc.h"
#include "CoulFunc.h"
#include <iostream>
#include <gsl/gsl_errno.h>

namespace Module_CrossSection{

    bool CrossSection::residualGamma_;
    bool CrossSection::residualNeutron_;
    bool CrossSection::residualProton_;
    bool CrossSection::residualAlpha_;
    bool CrossSection::calculateGammaCutoff_;
    std::vector<double> CrossSection::rateTemps_;
    std::vector<double> CrossSection::macsEnergies_;
    #endif
    bool Decayer::isCrossSection_;
    bool PreEqDecayer::isCrossSection_;
    double Decayer::maxL_;
    double PreEqDecayer::maxL_;
    double TransitionRateFunc::gammaCutoffEnergy_;
    ElementTable NuclearMass::elementTable_; 
    MassTable NuclearMass::massTable_;
    GDRTable GammaTransmissionFunc::gdrTable_;
    LevelsTable NuclearLevels::levelsTable_;
    int ParticleTransmissionFunc::alphaFormalism_;
    int ParticleTransmissionFunc::neutronFormalism_;
    int ParticleTransmissionFunc::protonFormalism_;
    int GammaTransmissionFunc::egdrType_;
    bool GammaTransmissionFunc::porterThomas_;
    bool ParticleTransmissionFunc::porterThomas_;

    void Initialize(); /**< Initialize default options*/

    void Go(int argc,char *argv[]); /**< Top level function to call from main*/
    void Run(int argc,char *argv[]); /**< Declaration of the main function of the CrossSection Module*/

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