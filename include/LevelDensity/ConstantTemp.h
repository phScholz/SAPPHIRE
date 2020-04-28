#pragma once

#include "LevelDensity/LevelDensityFormula.h"

class ConstantTemp : public LevelDensityFormula{
    public:
        ConstantTemp(int Z, int A, double J) : LevelDensityFormula(Z, A, J, 1){

        }
    
        double CalculateDensity(double E);
        double TotalLevelDensity(double E);
};