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
#include "Module_Decayer.h"
#include "SapphireInput.h"

//#include "Setup.cpp"

extern void Initialize();

/**
 * @brief Secondary function to print help. This is the first step to reconstruct Sapphire.cpp.
 */
void PrintHelp(){
  std::cout << std::endl;
  std::cout << " Supported modes:" << std::endl;
	std::cout << std::endl;
  std::cout << "\treaction      - Calculate reaction cross sections and/or rates." << std::endl;
  std::cout << "\tdecayer       - Calculate Monte-Carlo statistical decay." << std::endl;
  std::cout << std::endl;
  std::cout << "\told           - Instruction for the old Sapphire code." << std::endl;
	std::cout << std::endl;
  std::cout << "\thelp          - Show this help message." << std::endl;
  std::cout << "\ttemplate      - Print template input file." << std::endl;
	std::cout << std::endl;
  
} 

/**
 * @brief Check the first cmd parameter give and decide which module is responsible. Then run the respective module.
 * @param argc Number of cmd line parameters.
 * @param argv Array of cmd line parameters. 
 */
int main(int argc, char *argv[]) {
  /**
  * At the beginning of `main()` a header is printed and the cmd line arguments are checked.
  */
  std::cout << std::endl;
	std::cout << "Sapphire - A statistical nuclear reaction and decay code" << std::endl;
	std::cout << "********************************************************" << std::endl;
  std::cout << std::endl;

  /** If no argument is given at all, the PrintHelp() method is invoked and the program exits. */
  if (argc <=1 ) {
	  PrintHelp();
		return 0;
	}

  /** The Initializer function of Setup.cpp is called, to set some default values. */ 
  Initialize();

  /** The second cmd line parameter defines the module, which should be started, therefore it is stored in the std::string mode.*/
  std::string mode(argv[1]);
	
  /** 
  * Now it'll be checked if "mode" matches any of the available keywords: reaction, decayer, old, template, or help.
  * Depending on the matching, the respective Go() methods are called from Module_CrossSection.cpp, Module_Decayer.cpp, Module_OldSapphire.cpp, 
  * or the SapphireInput class.
  * If mode doesn't match any of the keywords, the PrintHelp() function is called.
  */
  if (mode == "reaction"){
    Module_CrossSection::Go(argc,argv);
  }
  else if (mode == "decayer"){
    Module_Decayer::Go(argc,argv);
  }
  else if (mode == "old"){
    Module_OldSapphire::Go(argc, argv);
  }
  else if (mode == "template"){
    SapphireInput *input = new SapphireInput();
    input->Go(argc,argv);
    delete input;
  }
  else if (mode == "help")
    PrintHelp();
   
  else
    PrintHelp();
  
  return 0;
}
