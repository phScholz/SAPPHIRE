#pragma once

#include "LevelDensity.h"
#include "NuclearMass.h"
#include <iostream>
#include <stdlib.h>
#include <math.h>

class RauscherLevelDensity : public LevelDensity {
 public:
  RauscherLevelDensity(int Z, int A, double J) : LevelDensity(Z,A,J) {
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


