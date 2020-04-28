#pragma once

#include "LevelDensity/LevelDensityFormula.h"

class GilbertCameron : public LevelDensityFormula{
    public:
        GilbertCameron(int Z, int A, double J) : LevelDensityFormula(Z, A, J, 1){

        }
    
        double CalculateDensity(double E);
        double TotalLevelDensity(double E);
};