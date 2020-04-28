/**
 * @file LevelDensityFormula.h
 * @brief Contains the LevelDensityFormula class
 * @date 2020-04-27
 */

#pragma once
#include "LevelDensity.h"

class LevelDensityFormula : public LevelDensity{
    public:
        LevelDensityFormula(int Z, int A, double J, int parity) : LevelDensity(Z, A, J, parity){
            SetTables(false); 
            criticalU_=2.5+150./A;
        }

    friend class KopeckyUhlGSF;
    friend class TransitionRateFunc;
    
    protected:
        virtual void CalcBackShift() = 0;
        virtual double CalcDensityParam(double) = 0;
        virtual double CalcNuclearTemp(double) = 0;
        void CalcConstantTempTerms();

    protected:
        double CalculateDensity(double E);
        double TotalLevelDensity(double E);
        
        

    protected:
        double backshift_;
        double criticalU_;
        double constAngTerm_;
        double nuclearTemp_;
        double e0_;
        static constexpr double zeta_ = 1.0;
        static constexpr double r0_ = 1.25;
};