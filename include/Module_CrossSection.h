/**
 * @file Module_CrossSection.h
 * @brief Reimplementation of the cross section calculations routine of Sapphire as a module.
 * @date 2020
 * @author Philipp Scholz, <pscholz@outlook.de>
 * 
 */

#ifndef MODULE_CROSSSECTION_H
#define MODULE_CROSSSECTION_H
#include <vector>
#include <map>
#include <string>

namespace Module_CrossSection{
    void Go(int argc,char *argv[]); /**< Top level function to call from main*/
    void Run(); /**< Declaration of the main function of the CrossSection Module*/
}