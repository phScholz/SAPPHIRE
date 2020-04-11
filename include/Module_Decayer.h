/**
 * @file Module_Decayer.h
 * @brief Reimplementation of the decay routine of Sapphire as a module.
 * @date 2020
 * @author Philipp Scholz, <pscholz@outlook.de>
 * 
 */

#ifndef MODULE_DECAYER_H
#define MODULE_DECAYER_H
#endif
#include <vector>
#include <map>
#include <string>

namespace Module_Decayer{
    void Go(int argc,char *argv[]); /**< Top level function to call from main*/
    void Run(int argc,char *argv[]); /**< Declaration of the main function of the Decayer Module*/

    /**
     * @brief Check wheter a string represents an actual file.
     * @param filename String whith the supposedly path to a file.
     * @return True if the file exists; False if it doesn't.
     */
    bool fexists(const char *filename);
}