/**
 * @file Sapphire.cpp
 * @brief Entry point for the Sapphire code. 
 * @author Philipp Scholz, <pscholz@outlook.de>
 * @date 2020
 * 
 * This version of Sapphire.cpp changed alot.
 * Now, for every task of Sapphire there should be a unique module, which should make the code more readable and understandable.
 * The project can be extended easily by simply adding new modules for new tasks.
 * Basically the whole content of the old and original Sapphire.cpp has been moved to Module_OldSapphire.cpp.
 */
#include <iostream>

/**
 * Now there are new includes which differ from the original version
 */
#include "Module_CrossSection.h"
#include "Module_OldSapphire.h"

/**
 * @brief Secondary function to print help. This is the first step to reconstruct Sapphire.cpp.
 */
void PrintHelp(){
  std::cout << std::endl;
	std::cout << "Statistical Analysis for Particle and Photoncapture and decay of HIgh energy REsonances" << std::endl;
	std::cout << "---------------------------------------------------------------------------------------" << std::endl;
	#ifndef MPI_BUILD
  std::cout << " MPI_BUILD is OFF" << std::endl;
  #endif 
  #ifdef MPI_BUILD
  std::cout << " MPI_BUILD is ON" << std::endl;
  #endif 
  std::cout << " Supported modes:" << std::endl;
	std::cout << std::endl;
  #ifndef MPI_BUILD
	std::cout << " reaction      - Calculate reaction cross sections and/or rates" << std::endl;
  #endif 
	std::cout << " decay         - Calculate Monte-Carlo statistical decay" << std::endl;
  std::cout << " old           - Instruction for the old Sapphire code." << std::endl;
	//std::cout << "Template     - Generate examplary configuration files" << std::endl;
	std::cout << " help           - Show this help message" << std::endl;
	std::cout << std::endl;
  
} 

int main(int argc, char *argv[]) {

  if (argc <=1 ) {
	  PrintHelp();
		return 0;
	}

  std::string mode(argv[1]);
	
  if (mode == "reaction")
		Module_CrossSection::Go(argc,argv);
  
  //else if (mode == "decay")

  else if (mode == "old")
    Module_oldSapphire::Go(argc, argv);
  
  else if (mode == "help")
    PrintHelp();
   
  else
    PrintHelp();
  
  return 0;
}