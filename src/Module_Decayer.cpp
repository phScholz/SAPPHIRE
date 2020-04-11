#include "Module_Decayer.h"
#include <iostream>
#include <fstream>

namespace Module_Decayer{
    void Go(int argc,char *argv[]){
       std::cout << "Not yet working."  << std::endl; 
    }
    
    void Run(int argc,char *argv[]){}

    bool fexists(const char *filename) {
        std::ifstream ifile(filename);
    return (bool)ifile;
    }
}