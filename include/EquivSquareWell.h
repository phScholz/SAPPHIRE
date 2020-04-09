#ifndef EQUIVSQUAREWELL_H
#define EQUIVSQUAREWELL_H

#include "CoulFunc.h"
#include "ParticleTransmissionFunc.h"

class TransmissionFunc;

class EquivSquareWell : public ParticleTransmissionFunc {
 public:
   EquivSquareWell(int z1, int m1, int z2, int m2, 
		   double jInitial, int piInitial,
		   double jFinal, int piFinal,
		   double spin, int parity, double maxL,
		   double totalWidthForCorrection,
		   double uncorrTotalWidthForCorrection,
		   double uncorrTotalWidthSqrdForCorrection,
		   TransmissionFunc* previous) :
    ParticleTransmissionFunc(z1,m1,z2,m2,jInitial,piInitial,
			     jFinal,piFinal,spin,parity,maxL,
			     totalWidthForCorrection,
			     uncorrTotalWidthForCorrection,
			     uncorrTotalWidthSqrdForCorrection,
			     previous) {
      coulFunc_ = new CoulFunc(z1_,z2_,redmass_,true);
   };    
   ~EquivSquareWell() {
      delete coulFunc_;
   };
   double CalcTransmission(double,int,double);
 private:
   double GetF();
   double GetV1();
   double GetV0();
   double GetR1();
   double GetR0();
 private:
  CoulFunc* coulFunc_;
};

#endif
