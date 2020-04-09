#ifndef MCCULLAGHGSF_H
#define MCCULLAGHGSF_H

#include "GammaTransmissionFunc.h"

class McCullaghGSF : public GammaTransmissionFunc {
 public:
  McCullaghGSF(int,int,double,int,double,int,double, 
	       double,double,double,TransmissionFunc*);
  double CalcStrengthFunction(double);
 private:
  double CalcEnergyDependentWidth(double,int);
};

#endif
