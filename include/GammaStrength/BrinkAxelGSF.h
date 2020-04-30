/** 
*   @file BrinkAxelGSF.h
*   @brief Header file for the Brink-Axel-GSF class
*   @date 2020
* 
*/
#pragma once
#ifndef BRINKAXELGSF_H
#define BRINKAXELGSF_H

#include "GammaTransmissionFunc.h"

/** @class BrinkAxelGSF
*   @brief Class for the Brink-Axel-GSF
*/
class BrinkAxelGSF : public GammaTransmissionFunc {
 public:

 /**
  * @brief Construct a BrinkAxelGSF as a child of a GammaTransmissionFunc.
  * @param z2 Charge/atomic number of the nucleus
  * @param m2 Mass number of the nucleus
  * @param jInitial Initial spin of the nucleus
  * @param piInitial Initial parity of the nucleus
  * @param jFinal Final spin of the nucleus
  * @param piFinal Final parity of the nucleus
  * @param maxL Maximum angular momentum
  * @param totalWidthForCorrection 
  * @param uncorrTotalWidthForCorrection
  * @param uncorrTotalWidthSqrdForCorrection
  * @param previous
  */
  BrinkAxelGSF(int z2, int m2, double jInitial, int piInitial,
			   double jFinal, int piFinal, double maxL, 
			   double totalWidthForCorrection,
			   double uncorrTotalWidthForCorrection,
			   double uncorrTotalWidthSqrdForCorrection,
			   TransmissionFunc* previous);

  /**
   * @brief Calculate the strength for a specific energy
   * @param energy Energy of a gamma-ray transition
   */
  double CalcStrengthFunction(double energy);

 private:

  double CalcEnergyDependentWidth(double,int);
};

#endif
