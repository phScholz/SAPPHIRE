#include "TransitionRateFunc.h"
#include "ParticleTransmissionFunc.h"
#include "GammaTransmissionFunc.h"
#include "RauscherLevelDensity.h"
#include "NuclearLevels.h"
#include <iostream>
#include <stdlib.h>
#include <omp.h>

TransitionRateFunc::TransitionRateFunc(int z1, int m1, int z2, int m2,
				       double jInitial, int piInitial, 
				       double jFinal, int piFinal,
				       double spin, int parity, double maxL,
				       double compoundE, double qValue,
				       double totalWidthForCorrection, 
				       double uncorrTotalWidthForCorrection, 
				       double uncorrTotalWidthSqrdForCorrection, 
				       TransitionRateFunc* previous,
				       bool isCrossSection) {
  levelDensity_ = new RauscherLevelDensity(z2,m2,jFinal);
  TransmissionFunc* previousTransmissionFunc = (previous) ? previous->transmissionFunc_ : NULL;
  if(m1!=0) {
    transmissionFunc_ =  
      ParticleTransmissionFunc::CreateParticleTransmissionFunc(z1,m1,z2,m2,
							       jInitial,piInitial,
							       jFinal,piFinal,
							       spin,parity,maxL,totalWidthForCorrection,
							       uncorrTotalWidthForCorrection,
							       uncorrTotalWidthSqrdForCorrection,
							       previousTransmissionFunc);
  } else {
    transmissionFunc_ = 
      GammaTransmissionFunc::CreateGammaTransmissionFunc(z2,m2,jInitial,piInitial,
							 jFinal,piFinal,maxL,
							 levelDensity_,totalWidthForCorrection,
							 uncorrTotalWidthForCorrection,
							 uncorrTotalWidthSqrdForCorrection,
							 previousTransmissionFunc,compoundE);
  }
  if(!transmissionFunc_->IsValid()) {
    std::cout << "An problem occurred creating transmission function. Aborting."
	      << std::endl;
    delete transmissionFunc_;
    exit(1);
  }
  double sum=0.;
  double exclusiveSum=0.;
  std::vector<Level> knownLevels = NuclearLevels::FindLevels(z2,m2);
  double highestKnownLevel = (knownLevels.size()) ? knownLevels[knownLevels.size()-1].energy_ : 0.;
  double lowEnergy=0.001; 
  double highEnergy = compoundE+qValue-highestKnownLevel;
  double dE  = (!isCrossSection) ? 0.01 : 0.05;
  double exclusiveLowEnergy = lowEnergy;
  if(isCrossSection) {
    double highestBoundEnergy=1000.;
    if(m1==0&&gammaCutoffEnergy_>0.) highestBoundEnergy = gammaCutoffEnergy_;
    else if(!NuclearMass::HighestBoundEnergy(z2,m2,highestBoundEnergy)) {
      //std::cout << "Could not calculate highest bound energy.  Aborting." << std::endl;
      //exit(1);
    }
    if(compoundE+qValue-highestBoundEnergy>exclusiveLowEnergy) 
      exclusiveLowEnergy = compoundE+qValue-highestBoundEnergy;
  }
  int numSteps = (highEnergy>lowEnergy) ? 50 : -50;//int((highEnergy-lowEnergy)/dE);
  if(numSteps>0) {
    if(numSteps%2!=0) {
      numSteps+=1;
    }
    dE = (highEnergy-lowEnergy)/double(numSteps);
  }
  double evenSum=0.;
  double oddSum=0.;
  double firstTerm=0.;
  double lastTerm=0.;

#pragma omp parallel for if(isCrossSection) reduction(+:sum,exclusiveSum,oddSum,evenSum)
  for(int i=0;i<=numSteps;++i) {
    double ep = lowEnergy+dE*i;
    double rate = CalcLevelDensity(compoundE+qValue-ep)*
      CalcTransmissionFunc(ep);
    if(i==0) {
      firstTerm = rate; 
    } else if(i==numSteps) {
      lastTerm = rate;
    } else {
      if(i%2!=0) {
	oddSum+=rate;
      } else {
	evenSum+=rate;
      }
    }
    if(ep<highEnergy) sum += rate*dE;
    if(ep>=exclusiveLowEnergy) exclusiveSum+=rate*dE;
    if(!isCrossSection&&ep<highEnergy){
      function_.push_back(XYPair(ep,rate*dE));
      cumulativeSum_.push_back(XYPair(ep+0.5*dE,sum));
    }
  }
  double continuousSum = dE/3.*(firstTerm + 2.*evenSum + 4.*oddSum + lastTerm);
  double discreteSum=0.;
  double groundStateTransmission = 0.;


#pragma omp parallel for if(isCrossSection) reduction(+:sum,exclusiveSum,discreteSum)
 for(int i=knownLevels.size()-1;i>=0;--i) {
    if(knownLevels[i].J_==jFinal&&knownLevels[i].Pi_==piFinal) {
      double ep = compoundE+qValue-knownLevels[i].energy_;
      if(ep<1e-10) continue;
      double rate = CalcTransmissionFunc(ep);
      discreteSum+=rate;
      sum+=rate;
      exclusiveSum+=rate;
      if(!isCrossSection){
	      function_.push_back(XYPair(ep,rate));
	      cumulativeSum_.push_back(XYPair(ep,sum));     
      }
      if(knownLevels[i].energy_==0.) groundStateTransmission = rate;
   }
  }
  integral_=continuousSum+discreteSum;
  groundStateTransmission_ = groundStateTransmission;
  if(!isCrossSection)
    for(int i = 0; i< cumulativeSum_.size(); i++) 
      if(integral_>0.) cumulativeSum_[i].Y_*=sum/integral_;
  exclusiveBranching_ = (sum>0.) ? exclusiveSum/sum : 0.;
}
