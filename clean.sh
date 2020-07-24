#!/bin/bash
##############################################
## @file	clean.sh
## @brief	Shell script to clean up all files which were created during a build routine
## @author	Philipp Scholz <pscholz@outlook.de>
## @date	2020-04-09
##############################################
# first, simply remove the build directory
rm -r build/ -f -v

# second, remove all the output in ./doc
rm -r ./doc/* -f -v 

# remove Makefiles and stuff
rm Makefile -v 
rm */Makefile -v
rm */CMakeCache.txt -v
rm CMakeCache.txt -v
rm coul/Makefile -v 
rm coul/cmake_install.cmake -v
rm coul/src/Makefile -v 
rm coul/src/cmake_install.cmake -v
rm coul/src/libcoul.a -v 
rm -r generated/ -v 
rm src/Makefile -v 
rm src/cmake_install.cmake -v
 
