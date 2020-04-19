#ifndef KOPECKYUHL_H
#define KOPECKYUHL_H
#pragma once
#include "GammaTransmissionFunc.h"

class LevelDensity;

class KopeckyUhlGSF : public GammaTransmissionFunc {
 public:
  KopeckyUhlGSF(int,int,double,int,double,int,double, 
		double,double,double,TransmissionFunc*,
		LevelDensity*,double);
  double CalcStrengthFunction(double);
 private:
  double CalcEnergyDependentWidth(double,int);
  double CalcKUTempSqrd(double);
 private:
  LevelDensity* levelDensity_;
  double compoundE_;
  double levelDenParam_;
};

#endif
