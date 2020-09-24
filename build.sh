#!/bin/bash
######################################################
## @file build.sh
## @author Philipp Scholz <pscholz@outlook.de>
## @date 2020-04-09
## @brief Script to build Sapphire
######################################################
#
echo ""
echo "##### SAPPHIRE #####"

#check if the build and output directory exists or not
#if not, then create the build directory
[ ! -d /build ] && mkdir -p ./build
[ ! -d /output ] && mkdir -p ./output

#go into the build directory
cd ./build

if [ -z $1 ]; then 
    buildtype="Debug";
else
    buildtype=$1;
fi


#run cmake and if sucessfull run make install 
cmake -DCMAKE_BUILD_TYPE=${buildtype} .. && make -j8 && make install

#run unit tests
if [ ${buildtype} != "Regulus" ]; then 
    echo ""; echo "##### RUNNING UNIT TESTS #####"; echo  ""; ./test_sapphire;
else
    echo ""; echo "##### UNIT TESTS ARE DISABLED #####"; echo "";
fi