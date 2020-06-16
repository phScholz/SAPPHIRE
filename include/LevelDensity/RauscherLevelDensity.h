/**
 * @file RauscherLevelDensity.h
 * @brief Contains RauscherLevelDensity as a child of LevelDensity.
 * @date 2020-04-22
 */

#pragma once

#include "LevelDensityFormula.h"
#include "Databases/NuclearMass.h"
#include <iostream>
#include <stdlib.h>
#include <math.h>

/**
 * @brief Class for the level density model by T. Rauscher. Child of LevelDensity.
 */
class RauscherLevelDensity : public LevelDensityFormula{
 public:
  RauscherLevelDensity(int Z, int A, double J) : LevelDensityFormula(Z,A,J,1) {
    SetTables(false);
    aTilda_ = alpha_*A+beta_*pow(A,0.666666667);
    microEnergyCorr_=0.;
    
    if(!NuclearMass::MicroEnergyCorr(Z,A,microEnergyCorr_)){
      //std::cout << "Microscopic energy correction not found.  Aborting." 
      //          << std::endl;
      //exit(1);
    }
    CalcBackShift();
    CalcConstantTempTerms();
    
  };

  RauscherLevelDensity(int Z, int A, double J, int P) : LevelDensityFormula(Z,A,J,P) {
    SetTables(false);
    aTilda_ = alpha_*A+beta_*pow(A,0.666666667);
    microEnergyCorr_=0.;
    
    if(!NuclearMass::MicroEnergyCorr(Z,A,microEnergyCorr_)){
      //std::cout << "Microscopic energy correction not found.  Aborting." 
      //          << std::endl;
      //exit(1);
    }
    CalcBackShift();
    CalcConstantTempTerms();
  };
  
  ~RauscherLevelDensity() {};
  void CalcBackShift();
  double CalcDensityParam(double);
  double CalcNuclearTemp(double);
  
 private:  
  static constexpr double alpha_ = 0.1337;
  static constexpr double beta_ = -0.06571;
  static constexpr double gamma_ = 0.04884;
  double aTilda_;
  double microEnergyCorr_;
};


