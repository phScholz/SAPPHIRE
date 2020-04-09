/*! \file BrinkAxelGSF.h
*   \brief Header file for the Brink-Axel-GSF class
* 
*/
#ifndef BRINKAXELGSF_H
#define BRINKAXELGSF_H

#include "GammaTransmissionFunc.h"

/*! \class BrinkAxelGSF
*   Class for the Brink-Axel-GSF
*/
class BrinkAxelGSF : public GammaTransmissionFunc {
 public:
  BrinkAxelGSF(int,int,double,int,double,int,double, 
	       double,double,double,TransmissionFunc*);
  double CalcStrengthFunction(double);
 private:
  double CalcEnergyDependentWidth(double,int);
};

#endif
