#include "GammaStrength/D1MQRPA.h"

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdexcept>

extern std::string sourceDirectory();

void D1MQRPA::ReadFile(){
    if(!FindE1())
        e1_=ReadE1();

    if(!FindM1())
        m1_=ReadM1();

    if(e1_){
        if(verbose_)PrintE1();
    }

    if(m1_){
        if(verbose_) PrintM1();
    }
}

bool D1MQRPA::FindE1(){
  QRPA_E1_Table::const_iterator it = e1Table.find(MassKey(z2_,m2_));

  if(it!=e1Table.end())
    { 
        e1Rows=it->second; 
        //std::cout << "DensityTable Entry found!" << std::endl;
        return true;
    }
  else
    {
        //std::cout << "Not found!" << std::endl;
        return false;
    }  
}

bool D1MQRPA::FindM1(){
    QRPA_M1_Table::const_iterator it = m1Table.find(MassKey(z2_,m2_));

  if(it!=m1Table.end())
    { 
        m1Rows=it->second; 
        //std::cout << "DensityTable Entry found!" << std::endl;
        return true;
    }
  else
    {
        //std::cout << "Not found!" << std::endl;
        return false;
    }  
}

bool D1MQRPA::ReadE1(){
    //construct filePaths
    int leading = 3;
    std::string fileE1 = sourceDirectory() + "/tables/gamma/d1m/z" + std::to_string(z2_*0.000001).substr(8-leading) + "_e1";
    if(verbose_) std::cout << "Reading E1 file ... " << fileE1 << std::endl;
    //construct search strings
    std::string massE1 = std::to_string(m2_) + " E1 PSF";
    if(verbose_) std::cout << "Search string ... " << massE1 << std::endl;

    //open file
    std::ifstream in(fileE1.c_str());

    //trying to read file
    if(!in) {
        throw std::runtime_error("Cannot read file "+fileE1);    
    }

    std::string line;

    //As long as we haven't found the right line to start read-in lines
    while( (line.find(massE1) == std::string::npos) && (!in.eof()) ){
        std::getline(in,line);
        //if(verbose_) std::cout << line.find(massE1) << std::endl;
    }

    if(in.eof()){
        throw std::runtime_error("Couldn't find mass "+ std::to_string(m2_)+ " in " +fileE1);
        return false;
    }
    if(verbose_) std::cout << "Found!!!" << std::endl;
    //skip one line
    for(int i = 0;i<1;++i) std::getline(in,line);    
    
    std::vector<double> dummyVector(10,0);

    for(int i=0; i<300; ++i){
        std::getline(in,line);
        std::istringstream lineStream(line);
        lineStream  >> dummyRowE1.energy;
        for(int i = 0; i < 10; i++){
            lineStream >> dummyVector.at(i);
        }
        dummyRowE1.strength = dummyVector;

        e1Rows.push_back(dummyRowE1);
    }
    in.close();
    e1Table[MassKey(z2_,m2_)] = e1Rows;
    return true;
}

bool D1MQRPA::ReadM1(){
    //construct filePaths
    int leading = 3;
    std::string fileM1 = sourceDirectory() + "/tables/gamma/d1m/z" + std::to_string(z2_*0.000001).substr(8-leading) + "_m1";
    if(verbose_) std::cout << "Reading M1 file ... " << fileM1 << std::endl;
    //construct search strings
    std::string massM1 = std::to_string(m2_) + " M1 PSF";

    //open file
    std::ifstream in(fileM1.c_str());

    //trying to read file
    if(!in) {
        throw std::runtime_error("Cannot read file "+fileM1);    
    }

    std::string line;

    //As long as we haven't found the right line to start read-in lines
    while( (line.find(massM1) == std::string::npos) && (!in.eof() )){
        std::getline(in,line);
    }

    if(in.eof()){
        throw std::runtime_error("Couldn't find mass "+ std::to_string(m2_)+ " in " +fileM1);
        return false;
    }
    if(verbose_) std::cout << "Found!!!" << std::endl;

    //skip one line
    for(int i = 0;i<1;++i) std::getline(in,line);    
    
    std::vector<double> dummyVector(10,0);

    for(int i=0; i<300; ++i){
        std::getline(in,line);
        std::istringstream lineStream(line);
        lineStream  >> dummyRowM1.energy;
        for(int i = 0; i < 2; i++){
            lineStream >> dummyVector.at(i);
        }
        dummyRowM1.strength = dummyVector;

        m1Rows.push_back(dummyRowM1);
    }
    m1Table[MassKey(z2_,m2_)] = m1Rows  ;
    in.close();
    return true;
}

void D1MQRPA::PrintE1(){
    if(verbose_) std::cout << std::endl << "Printing E1 strength " << std::endl;
    std::cout.precision(4);
    for(std::vector<QRPAE1row>::iterator it = e1Rows.begin(); it !=e1Rows.end(); ++it){
        std::cout << std::scientific << it->energy << "\t" << it->strength.at(0) << std::endl;
    }
}

void D1MQRPA::PrintM1(){
    if(verbose_) std::cout << std::endl << "Printing M1 strength " << std::endl;
    std::cout.precision(4);
    for(std::vector<QRPAM1row>::iterator it = m1Rows.begin(); it !=m1Rows.end(); ++it){
        std::cout << std::scientific <<it->energy << "\t" <<it->strength.at(0) << std::endl;
    }
}

double D1MQRPA::CalcStrengthFunction(double energy){

    //opposite parities && dL = 1
    if(piFinal_==-1*piInitial_ && abs(jFinal_-jInitial_) == 1){
        return CalcE1Strength(energy);
    }

    // same parities && dL <=1 && (one of the two spins need to be not zero)
    if((piFinal_==piInitial_) && (abs(jFinal_-jInitial_) <= 1) && (jFinal_!=0 || jInitial_!=0) ){
        return CalcM1Strength(energy);
    }
}

double D1MQRPA::CalcE1Strength(double energy){

    if(e1Rows.size()>0){
        //Looking for the energy in the e1Rows
        double lowerx=0;
        double lowery=0;
        double upperx=0;
        double uppery=0;

        double strength=0;
    

        for(auto it = e1Rows.begin(); it != e1Rows.end(); ++it){
            if (it->energy == energy){
                strength = it->strength.at(0);
                return strength;
            }

            //Linear Interpolation
            if(it->energy > lowerx && it->energy < energy){
                lowerx = it->energy;
                lowery = it->strength.at(0);

            }
            
            if(it->energy > energy && upperx < energy){
                upperx = it->energy;
                uppery = it->strength.at(0);
            }
        }

        //If the energy doenst exist in e1Rows we need to interpolate
        strength =(energy-lowerx)*(uppery-lowery)/(upperx-lowerx)+lowery;
        
        return strength+f0_*exEnergy_/(1+exp(energy-e0_));
    }
}

double D1MQRPA::CalcM1Strength(double energy){

    if(m1Rows.size()>0){
        //Looking for the energy in the e1Rows
        double lowerx=0;
        double lowery=0;
        double upperx=0;
        double uppery=0;

        double strength=0;
    

        for(auto it = m1Rows.begin(); it != m1Rows.end(); ++it){
            if (it->energy == energy){
                strength = it->strength.at(0);
                return strength;
            }

            //Linear Interpolation
            if(it->energy > lowerx && it->energy < energy){
                lowerx = it->energy;
                lowery = it->strength.at(0);

            }
            
            if(it->energy > energy && upperx < energy){
                upperx = it->energy;
                uppery = it->strength.at(0);
            }
        }

        //If the energy doenst exist in e1Rows we need to interpolate
        strength =(energy-lowerx)*(uppery-lowery)/(upperx-lowerx)+lowery;
        
        return strength+C_*exp(-1.0*eta_*energy);
    }
}

void D1MQRPA::SetUpperLimit(){
    f0_ = 5.0e-10;
    e0_ = 5.0;
    eta_ = 0.8;
    C_ = 3.0e-8;
}

void D1MQRPA::SetLowerLimit(){
    f0_ = 1.0e-10;
    e0_ = 3.0;
    eta_ = 0.8;
    C_ = 1.0e-8;
}