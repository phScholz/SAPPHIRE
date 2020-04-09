#include "McCullaghGSF.h"

McCullaghGSF::McCullaghGSF(int z2, int m2, double jInitial, int piInitial,
			   double jFinal, int piFinal, double maxL, 
			   double totalWidthForCorrection,
			   double uncorrTotalWidthForCorrection,
			   double uncorrTotalWidthSqrdForCorrection,
			   TransmissionFunc* previous) :
  GammaTransmissionFunc(z2,m2,jInitial,piInitial,jFinal,piFinal,maxL,
			totalWidthForCorrection,uncorrTotalWidthForCorrection,
			uncorrTotalWidthSqrdForCorrection,previous) {
}

double McCullaghGSF::CalcStrengthFunction(double energy) {
  double strength = 0.;
  for(int i = 0;i<2;i++) {
    strength+=gdrParameters_.kSigmaGamma_[i]*CalcEnergyDependentWidth(energy,i)/
      (pow(pow(energy,2.)-pow(gdrParameters_.E_[i],2.),2.)+
       pow(CalcEnergyDependentWidth(energy,i)*energy,2.));
  }
  if(maxL_==0.||maxL_==1.) strength*=energy;
  else strength/=energy;
  return strength;
}

double McCullaghGSF::CalcEnergyDependentWidth(double energy, int which) {
  double width = gdrParameters_.W_[which]*sqrt(energy/gdrParameters_.E_[which]);
  return width;
}
