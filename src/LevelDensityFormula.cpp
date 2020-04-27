/**
 * @file LevelDensity.cpp
 * @brief Logic of the LevelDensityFormula class
 * @date 2020-04-27
 */

#include "LevelDensityFormula.h"
#include "Constants.h"
#include <math.h>

double LevelDensityFormula::CalculateDensity(double E){
    double density;
    
    if(E-backshift_>criticalU_) {
      double u = E-backshift_;
      double a = CalcDensityParam(u);
      double q = sqrt(2./5.*r0_*r0_*uconv/hbarc/hbarc*zeta_)*
        pow(A_,0.83333333333)*pow(u/a,0.25);

      double angTerm = (2.*J_+1.)/4./q/q*exp(-(J_+0.5)*(J_+0.5)/2./q/q);
      //double angTerm = (2.*J_+1.)/4./q/q*exp(-J_*(J_+1.)/2./q/q);
      double totalDensity = 1/12./sqrt(2.)/q/pow(a,0.25)*exp(2.*sqrt(a*u))/
        pow(u,1.25);

      density = angTerm*totalDensity;
    } else {
      double totalDensity = 1./nuclearTemp_*exp((E-e0_)/nuclearTemp_);

    density = constAngTerm_*totalDensity;
    }

    return density;
}

double LevelDensityFormula::TotalLevelDensity(double E) {
  double density;
  if(E-backshift_>criticalU_) {
    double u = E-backshift_;
    double a = CalcDensityParam(u);
    double q = sqrt(2./5.*r0_*r0_*uconv/hbarc/hbarc*zeta_)*
      pow(A_,0.83333333333)*pow(u/a,0.25);
    density  = 1/12./sqrt(2.)/q/pow(a,0.25)*exp(2.*sqrt(a*u))/
      pow(u,1.25);
  } else {
    density = 1./nuclearTemp_*exp((E-e0_)/nuclearTemp_);
  }
  return density;
}

void LevelDensityFormula::CalcConstantTempTerms() {
  double a = CalcDensityParam(criticalU_);
  double q = sqrt(2./5.*r0_*r0_*uconv/hbarc/hbarc*zeta_)*
    pow(A_,0.83333333333)*pow(criticalU_/a,0.25);  
 
  constAngTerm_ = (2.*J_+1.)/4./q/q*exp(-(J_+0.5)*(J_+0.5)/2./q/q);
  //constAngTerm_ = (2.*J_+1.)/4./q/q*exp(-J_*(J_+1.)/2./q/q);
  nuclearTemp_=CalcNuclearTemp(criticalU_);

  double totalDensity = 1/12./sqrt(2.)/q/pow(a,0.25)*
    exp(2.*sqrt(a*criticalU_))/pow(criticalU_,1.25);
  e0_ = criticalU_+backshift_-nuclearTemp_*log(nuclearTemp_*totalDensity);
}