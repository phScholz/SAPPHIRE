/**
 * @file LevelDensityTable.h
 * @date 2020-04-27
 * @brief Declaration and Definition of the LevelDensityTable class
 */

#pragma once
#include "LevelDensity.h"
#include <string>
#include <vector>
#include <iostream>

class LevelDensityTable : public LevelDensity{
    
    protected:
        virtual void GetFileName() = 0;
        virtual void ReadFile() = 0;
        virtual void FillVector() = 0;
        virtual double TotalLevelDensity(double E) = 0;
        virtual double CalculateDensity(double E) = 0;
        void SetFilename(std::string x){filename=x;}

    public:
        LevelDensityTable(int Z, int A, double J, int parity): LevelDensity(Z,A,J,parity){

        }

    protected:
        std::vector< std::pair<double, double> > DensityVector; /**< Densities will be read to this vector*/
        std::string filename; /**< Filename of the data table*/
};
