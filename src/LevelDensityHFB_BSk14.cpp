/**
 * @file LevelDensityHFB_BSk14.cpp
 * @brief Contains the logic for the HFB-BSk14 datatables
 * @date 2020-04-27
 */
#include "LevelDensity/LevelDensityHFB_BSk14.h"
#include "NuclearMass.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <boost/math/interpolators/barycentric_rational.hpp>

extern std::string sourceDirectory();

void LevelDensityHFB_BSk14::GetFileName(){
    if(verbose_) std::cout << "Getting FileName ... "  << std::endl;
    std::string element = NuclearMass::FindElement(Z_);
    SetFilename(sourceDirectory() + "/tables/densities/HFB-2008/" + element + ".tab");
}

void LevelDensityHFB_BSk14::ReadFile(){
    if(verbose_) std::cout << "Reading File " << filename << std::endl;
    //define search strings
    std::string mass = std::to_string(A_)+":";
    if(verbose_) std::cout << "Mass string: " << mass << std::endl;
    std::string parity;
    if(parity_>0) parity = "Positive-Parity";
    else parity = "Negative-parity";

    if(verbose_) std::cout << "Parity string: " << parity << std::endl;

    //open file
    std::ifstream in(filename.c_str());
    
    //trying to read file
    if(!in) {
        throw std::runtime_error("Cannot read file "+filename);    
    }

    std::string line;

    //As long as we haven't found the right line to start read-in lines
    while( ((line.find(mass) == std::string::npos) || (line.find(parity) == std::string::npos)) && (!in.eof())){
        std::getline(in,line);
        if(verbose_) std::cout << line.find(mass) << " " << line.find(parity) << std::endl;
    }
    if(verbose_){
        if(line.find(mass) != std::string::npos){
            if(line.find(parity) != std::string::npos){
                std::cout << "Found" << std::endl;
                std::cout << line << std::endl;
            }
        }
    }
    

    if(in.eof()){
        throw std::runtime_error("Couldn't find mass "+ std::to_string(A_)+ " in " +filename);
    }
    
    //skip the first four lines
    for(int i = 0;i<4;++i) std::getline(in,line);    
    //create a dummy row which will then pushed back to the rows vector
    
    for(int i=0; i<56; ++i){
        std::getline(in,line);
        std::istringstream lineStream(line);
        lineStream  >> dummyRow.U
                    >> dummyRow.T
                    >> dummyRow.NCUMUL
                    >> dummyRow.RHOOBS 
                    >> dummyRow.RHOTOT 
                    >> dummyRow.J0
                    >> dummyRow.J1
                    >> dummyRow.J2
                    >> dummyRow.J3
                    >> dummyRow.J4
                    >> dummyRow.J5
                    >> dummyRow.J6
                    >> dummyRow.J7
                    >> dummyRow.J8
                    >> dummyRow.J9
                    >> dummyRow.J10
                    >> dummyRow.J11
                    >> dummyRow.J12
                    >> dummyRow.J13
                    >> dummyRow.J14
                    >> dummyRow.J15
                    >> dummyRow.J16
                    >> dummyRow.J17
                    >> dummyRow.J18
                    >> dummyRow.J19
                    >> dummyRow.J20
                    >> dummyRow.J21
                    >> dummyRow.J22
                    >> dummyRow.J23
                    >> dummyRow.J24
                    >> dummyRow.J25
                    >> dummyRow.J26
                    >> dummyRow.J27
                    >> dummyRow.J28
                    >> dummyRow.J29
                    >> dummyRow.J30
                    >> dummyRow.J30
                    >> dummyRow.J31
                    >> dummyRow.J32
                    >> dummyRow.J33
                    >> dummyRow.J34
                    >> dummyRow.J35
                    >> dummyRow.J36
                    >> dummyRow.J37
                    >> dummyRow.J38
                    >> dummyRow.J39
                    >> dummyRow.J40
                    >> dummyRow.J40
                    >> dummyRow.J41
                    >> dummyRow.J42
                    >> dummyRow.J43
                    >> dummyRow.J44
                    >> dummyRow.J45
                    >> dummyRow.J46
                    >> dummyRow.J47
                    >> dummyRow.J48
                    >> dummyRow.J49;

        rows.push_back(dummyRow);
    }
    in.close();
    densityTable[ChargeMassParityKey(Z_,A_,parity_)]=rows;    
}

void LevelDensityHFB_BSk14::PrintDensityTable(){
    for(auto it = densityTable.begin(); it!=densityTable.end(); ++it){
        std::cout << it->first.A_ << " " << it->first.Z_ << " " << it->first.P_ << std::endl;
    }
}

void LevelDensityHFB_BSk14::FillVector(){
    for (auto it = rows.begin(); it !=rows.end(); ++it){
        switch ((int) floor(J_)){
            case 0: DensityVector.push_back(std::pair<double,double>(it->U, it->J0)); break;
            case 1: DensityVector.push_back(std::pair<double,double>(it->U, it->J1)); break; 
            case 2: DensityVector.push_back(std::pair<double,double>(it->U, it->J2)); break;
            case 3: DensityVector.push_back(std::pair<double,double>(it->U, it->J3)); break;
            case 4: DensityVector.push_back(std::pair<double,double>(it->U, it->J4)); break;
            case 5: DensityVector.push_back(std::pair<double,double>(it->U, it->J5)); break;        
            case 6: DensityVector.push_back(std::pair<double,double>(it->U, it->J6)); break;
            case 7: DensityVector.push_back(std::pair<double,double>(it->U, it->J7)); break;
            case 8: DensityVector.push_back(std::pair<double,double>(it->U, it->J8)); break;
            case 9: DensityVector.push_back(std::pair<double,double>(it->U, it->J9)); break;
            case 10:DensityVector.push_back(std::pair<double,double>(it->U, it->J10)); break;
            case 11:DensityVector.push_back(std::pair<double,double>(it->U, it->J11)); break;
            case 12:DensityVector.push_back(std::pair<double,double>(it->U, it->J12)); break;
            case 13:DensityVector.push_back(std::pair<double,double>(it->U, it->J13)); break;
            case 14:DensityVector.push_back(std::pair<double,double>(it->U, it->J14)); break;
            case 15:DensityVector.push_back(std::pair<double,double>(it->U, it->J15)); break;
            case 16:DensityVector.push_back(std::pair<double,double>(it->U, it->J16)); break;
            case 17:DensityVector.push_back(std::pair<double,double>(it->U, it->J17)); break;
            case 18:DensityVector.push_back(std::pair<double,double>(it->U, it->J18)); break;
            case 19:DensityVector.push_back(std::pair<double,double>(it->U, it->J19)); break;
            case 20:DensityVector.push_back(std::pair<double,double>(it->U, it->J20)); break;
            case 21:DensityVector.push_back(std::pair<double,double>(it->U, it->J21)); break;
            case 22:DensityVector.push_back(std::pair<double,double>(it->U, it->J22)); break;
            case 23:DensityVector.push_back(std::pair<double,double>(it->U, it->J23)); break;
            case 24:DensityVector.push_back(std::pair<double,double>(it->U, it->J24)); break;
            case 25:DensityVector.push_back(std::pair<double,double>(it->U, it->J25)); break;
            case 26:DensityVector.push_back(std::pair<double,double>(it->U, it->J26)); break;
            case 27:DensityVector.push_back(std::pair<double,double>(it->U, it->J27)); break;
            case 28:DensityVector.push_back(std::pair<double,double>(it->U, it->J28)); break;
            case 29:DensityVector.push_back(std::pair<double,double>(it->U, it->J29)); break;
            case 30:DensityVector.push_back(std::pair<double,double>(it->U, it->J30)); break;
            case 31:DensityVector.push_back(std::pair<double,double>(it->U, it->J31)); break;
            case 32:DensityVector.push_back(std::pair<double,double>(it->U, it->J32)); break;
            case 33:DensityVector.push_back(std::pair<double,double>(it->U, it->J33)); break;
            case 34:DensityVector.push_back(std::pair<double,double>(it->U, it->J34)); break;
            case 35:DensityVector.push_back(std::pair<double,double>(it->U, it->J35)); break;
            case 36:DensityVector.push_back(std::pair<double,double>(it->U, it->J36)); break;
            case 37:DensityVector.push_back(std::pair<double,double>(it->U, it->J37)); break;
            case 38:DensityVector.push_back(std::pair<double,double>(it->U, it->J38)); break;
            case 39:DensityVector.push_back(std::pair<double,double>(it->U, it->J39)); break;
            case 40:DensityVector.push_back(std::pair<double,double>(it->U, it->J40)); break;
            case 41:DensityVector.push_back(std::pair<double,double>(it->U, it->J41)); break;
            case 42:DensityVector.push_back(std::pair<double,double>(it->U, it->J42)); break;
            case 43:DensityVector.push_back(std::pair<double,double>(it->U, it->J43)); break;
            case 44:DensityVector.push_back(std::pair<double,double>(it->U, it->J44)); break;
            case 45:DensityVector.push_back(std::pair<double,double>(it->U, it->J45)); break;
            case 46:DensityVector.push_back(std::pair<double,double>(it->U, it->J46)); break;
            case 47:DensityVector.push_back(std::pair<double,double>(it->U, it->J47)); break;
            case 48:DensityVector.push_back(std::pair<double,double>(it->U, it->J48)); break;
            case 49:DensityVector.push_back(std::pair<double,double>(it->U, it->J49)); break;

            default:
                throw std::runtime_error("Spin is to large for NLD tables "+ std::to_string(J_)+ "!!!");
        }
    }

}

void LevelDensityHFB_BSk14::PrintRows(){
    if(verbose_) std::cout << "Printing Rows ... "  << std::endl;
    std::cout << std::endl;

    std::cout   
            << "U" << "        "
            << "T" << "        "
            << "NCUMUL" << "   "
            << "RHOOBS" << "   "
            << "RHOTOT" << "   "
            << "J0" << "       "
            << "J1" << "       "
            << "J2" << "       "
            << "J3" << "       "
            << "J4" << "       "
            << "J5" << "       "
            << "J6" << "       "
            << "J7" << "       "
            << "J8" << "       "
            << "J9" << "       "
            << "J10" << std::endl;
    
    for(std::vector<HFBTabRow>::iterator it = rows.begin(); it !=rows.end(); ++it){
        it->PrintRow();
    }
}

double LevelDensityHFB_BSk14::CalculateDensity(double E){
    if(DensityVector.size()>0){
        //Looking for the energy in the DensityVector
        double lowerx=0;
        double lowery=0;
        double upperx=0;
        double uppery=0;
        for(auto it = DensityVector.begin(); it != DensityVector.end(); ++it){
            if (it->first == E){
                return it->second;
            }
            if(it->first > lowerx && it->first < E){
                lowerx = it->first;
                lowery = it->second;

            }
            
            if(it->first > E && upperx < E){
                upperx = it->first;
                uppery = it->second;
            }
        }

        //If the energy doenst exist in the DensityVector we need to interpolate

        return (E-lowerx)*(uppery-lowery)/(upperx-lowerx)+lowery;
    }
    else
    {
        throw std::runtime_error("No entries in DensityVector !!!");
        return 0;
    }
}

double LevelDensityHFB_BSk14::TotalLevelDensity(double E){
    if(rows.size()>0){
        //Looking for the energy in the DensityVector
        double lowerx=0;
        double lowery=0;
        double upperx=0;
        double uppery=0;
        for(auto it = rows.begin(); it != rows.end(); ++it){
            if (it->U == E){
                return it->RHOTOT;
            }
            if(it->U > lowerx && it->U < E){
                lowerx = it->U;
                lowery = it->RHOTOT;

            }
            
            if(it->U > E && upperx < E){
                upperx = it->U;
                uppery = it->RHOTOT;
            }
        }

        //If the energy doenst exist in the DensityVector we need to interpolate

        return (E-lowerx)*(uppery-lowery)/(upperx-lowerx)+lowery;
    }
    else
    {
        throw std::runtime_error("No entries in DensityVector !!!");
        return 0;
    }
}

bool LevelDensityHFB_BSk14::FindDensities(int Z, int A, int P){    
  HFBTable::const_iterator it = densityTable.find(ChargeMassParityKey(Z,A,P));

  if(it!=densityTable.end())
    { 
        rows=it->second; 
        //std::cout << "DensityTable Entry found!" << std::endl;
        return true;
    }
  else
    {
        //std::cout << "Not found!" << std::endl;
        return false;
    }  
}


void HFBTabRow::PrintRow(){
        std::cout.precision(1);
        std::cout   << std::scientific
                    << U << "  "
                    << T << "  "
                    << NCUMUL << "  "
                    << RHOOBS << "  "
                    << RHOTOT << "  "
                    << J0 << "  "
                    << J1 << "  "
                    << J2 << "  "
                    << J3 << "  "
                    << J4 << "  "
                    << J5 << "  "
                    << J6 << "  "
                    << J7 << "  "
                    << J8 << "  "
                    << J9 << "  "
                    << J10 << std::endl;
}




