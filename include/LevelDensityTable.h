/**
 * @file LevelDensityTable.h
 * @date 2020-04-27
 * @brief Declaration and Definition of the LevelDensityTable class
 */

#pragma once
#include "LevelDensity.h"
#include <string>
#include <vector>
#include <tuple>

class LevelDensityTable : public LevelDensity{
    public:
        LevelDensityTable(int Z, int A, double J, int parity): LevelDensity(Z,A,J,parity){
            SetTables(true);
            GetFileName(Z_,A_);
            FillVector(J_,parity_);
        }
    
    protected:
    virtual void GetFileName(int Z, int A){};
    virtual void FillVector(double J, int parity){}; 

    protected:
        std::string filename; /**< Filename of the data table*/
};
