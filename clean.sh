#!/usr/bin/bash
##############################################
## @file	clean.sh
## @brief	Shell script to clean up all files which were created during a build routine
## @author	Philipp Scholz <pscholz@outlook.de>
## @date	2020-04-09
##############################################
# first, simply remove the build directory
rm -r build/

# second, remove all the output in ./doc
rm -r ./doc/*

# remove Makefiles and stuff
rm Makefile
rm coul/Makefile
rm coul/cmake_install.cmake
rm coul/src/Makefile
rm coul/src/cmake_install.cmake
rm coul/src/libcoul.a
rm -r generated/
rm src/Makefile
rm src/cmake_install.cmake
 
