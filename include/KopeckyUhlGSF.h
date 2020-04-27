/**
 * @file KopeckyUhlGSF.h
 * @brief Declaration for the class of the GLO
 * @date 2020-04-27
 */
#pragma once
#include "GammaTransmissionFunc.h"
#include "LevelDensityFormula.h"


class KopeckyUhlGSF : public GammaTransmissionFunc {
 public:
  KopeckyUhlGSF(int,int,double,int,double,int,double, 
		double,double,double,TransmissionFunc*,
		LevelDensityFormula*,double);
  double CalcStrengthFunction(double);
 private:
  double CalcEnergyDependentWidth(double,int);
  double CalcKUTempSqrd(double);
 private:
  LevelDensityFormula* levelDensity_;
  double compoundE_;
  double levelDenParam_;
};


