/**
 * @file Module_BreitWigner.h
 * @brief Contains methods to calculate HF and BW cross sections from the input of beta-delayed neutron experiments
 */

#pragma once
#include <iostream>

namespace Module_BreitWigner{
    /** 
    *   @brief Top level function to call from main for the reaction Module.
    *   @param argc Number of cmd line parameters of Sapphire.
    *   @param argv The array which contains the cmd line parameters of Sapphire.
    */
    void Go(int argc,char *argv[]); 

    /**
     * @brief Printing help information to std::cout
     */
    void PrintHelp();

    


}