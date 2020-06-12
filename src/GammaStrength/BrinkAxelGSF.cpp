#include "GammaStrength/BrinkAxelGSF.h"

double BrinkAxelGSF::CalcStrengthFunction(double energy) {
  double strength = 0.;
  
  if(energy == 0){
    return 0;
  }
  else{

    for(int i = 0;i<2;i++) {
      strength+=gdrParameters_.kSigmaGamma_[i]*CalcEnergyDependentWidth(energy,i)/
        (pow(pow(energy,2.)-pow(gdrParameters_.E_[i],2.),2.)+
         pow(CalcEnergyDependentWidth(energy,i)*energy,2.));
    }

    if(maxL_==0.|| maxL_==1.){
       strength*=energy;
    }
    else{
      strength/=energy;
    }

    return strength;
  }
}

double BrinkAxelGSF::CalcEnergyDependentWidth(double energy, int which) {
  double width = gdrParameters_.W_[which];
  return width;
}
