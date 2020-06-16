/**
 * @file Module_BetaNeutron.h
 * @brief Contains methods to calculate HF and BW cross sections from the input of beta-delayed neutron experiments
 */

#pragma once

namespace Module_BetaNeutron{
    /** 
    *   @brief Top level function to call from main for the reaction Module.
    *   @param argc Number of cmd line parameters of Sapphire.
    *   @param argv The array which contains the cmd line parameters of Sapphire.
    */
    void Go(int argc,char *argv[]); 

    /**
     * @brief Fit the model for the neutron-OMP to the given partial neutron widths.
     */
    void FitNeutronWidth(){};

    /**
     * @brief Fit the level density model to the given level spacing.
     */
    void FitLevelDensity(){};

    /**
     * @brief Fit gamma strength model to given partial gamma widths.
     */
    void FitGammaStrength(){};

    /**
     * @brief Calculate HF cross section for fitted parameters.
     */
    void CalculateHF(){};

    /**
     * @brief Calculate BW cross section for fitted parameters.
     */
    void CalculateBW(){};

}