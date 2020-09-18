#!/bin/bash
######################################################
## @file build.sh
## @author Philipp Scholz <pscholz@outlook.de>
## @date 2020-04-09
## @brief Script to build Sapphire
######################################################
#
echo ""
echo "##### BUILDING SAPPHIRE #####"
echo  ""

#check if the build and output directory exists or not
#if not, then create the build directory
[ ! -d /build ] && mkdir -p ./build
[ ! -d /output ] && mkdir -p ./output

#go into the build directory
cd ./build

#run cmake and if sucessfull run make install 
cmake -DCMAKE_BUILD_TYPE=$1 .. && make -j8 && make install

echo ""
echo "##### RUNNING UNIT TESTS #####"
echo  ""
#run unit tests
./test_sapphire 

