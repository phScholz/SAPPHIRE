/**
 * @file KopeckyUhlGSF.h
 * @brief Declaration for the class of the GLO
 * @date 2020-04-27
 */
#pragma once
#include "GammaTransmissionFunc.h"
#include "LevelDensity/RauscherLevelDensity.h"


class KopeckyUhlGSF : public GammaTransmissionFunc {
 public:
  KopeckyUhlGSF(int,int,double,int,double,int,double, 
		double,double,double,TransmissionFunc*, double);
  double CalcStrengthFunction(double);
 private:
  double CalcEnergyDependentWidth(double,int);
  double CalcKUTempSqrd(double);
 private:
  RauscherLevelDensity* levelDensity_;
  double compoundE_;
  double levelDenParam_;
};


