/**
 * @file Module_RandomScheme.h
 * @date 2020-04-25
 * @brief Module to create random level schemes and calculate features.
 * 
 */
#pragma once
#include "SapphireInput.h"

#include <iostream>
#include <fstream>
#include <string>

namespace Module_RandomScheme{
    void Go(int argc,char *argv[]); /**< Top level function to call from main*/
    void RunCreate(const SapphireInput & input); /**< Declaration of the main function of the Decayer Module*/
    void RunExtend(const SapphireInput & input); /**< Declaration of the main function of the Decayer Module*/
    void PrintHelp();
    void Create(int argc, char *argv[]);
    void Extend(int argc, char *argv[]);
    
    /**
     * @brief Check wheter a string represents an actual file.
     * @param filename String with the supposedly path to a file.
     * @return True if the file exists; False if it doesn't.
     */
    bool fileExists(const char *filename);
}

