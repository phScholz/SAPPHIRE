#include "Decayer.h"
#include "TransitionRateFunc.h"
#include "NuclearMass.h"
#include "NuclearLevels.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#ifndef MPI_BUILD
#include <omp.h>
#endif

#ifndef MPI_BUILD
extern unsigned int randomSeed[12];
#endif

Decayer::Decayer(int Z, int A, double jInitial, int piInitial, double energy,
		 double totalWidthForCorrection,
		 double uncorrTotalWidthForCorrection,
		 double uncorrTotalWidthSqrdForCorrection,
		 Decayer* widthCorrectedDecayer) :
  Z_(Z), A_(A), piInitial_(piInitial), jInitial_(jInitial), energy_(energy),
  totalWidthForCorrection_(totalWidthForCorrection), 
  uncorrTotalWidthForCorrection_(uncorrTotalWidthForCorrection),
  uncorrTotalWidthSqrdForCorrection_(uncorrTotalWidthSqrdForCorrection),
  widthCorrectedDecayer_(widthCorrectedDecayer) {

  neutronEntrance_=0.;
  gammaEntrance_=0.;
  protonEntrance_=0.;
  alphaEntrance_=0.;
  neutronTotalWidth_=0.;
  alphaTotalWidth_=0.;
  gammaTotalWidth_=0.;
  protonTotalWidth_=0.;

  //Check if state is a known bound state
  double qValueProton,qValueNeutron,qValueAlpha;
  if(!NuclearMass::QValue(Z,A,Z,A-1,qValueNeutron)) {
    std::cout << "Unknown masses requested.  Aborting." << std::endl;
    exit(1);
  }
  if(!NuclearMass::QValue(Z,A,Z-1,A-1,qValueProton)) {
    std::cout << "Unknown masses requested.  Aborting." << std::endl;
    exit(1);
  }
  if(!NuclearMass::QValue(Z,A,Z-2,A-4,qValueAlpha)) {
    std::cout << "Unknown masses requested.  Aborting." << std::endl;
    exit(1);
  }
  std::vector<Level> knownLevels = NuclearLevels::FindLevels(Z,A);
  std::vector<Level>::const_iterator foundLevel = knownLevels.end();
  if(qValueAlpha+energy<=0.&&qValueProton+energy<=0.&&qValueNeutron+energy<=0.) {
    for(std::vector<Level>::const_iterator it = knownLevels.begin();
	it<knownLevels.end();++it) {
      if((it->energy_-1.e-10<=energy&&energy<=it->energy_+1.e-10)&&
	 jInitial==it->J_&&piInitial==it->Pi_) {
	foundLevel = it;
	break;
      }
    }  
  }
  totalIntegral_=0.;    
  totalIntegralSqrd_=0.;  
  bool knownCDFBuilt = false;
  if(foundLevel!=knownLevels.end()&&foundLevel->gammas_.size()>0) {
    //if level is known, and has a known decay scheme use it
     knownCDFBuilt = BuildKnownCDF(int(foundLevel-knownLevels.begin()),knownLevels);
  } 
  if(!knownCDFBuilt) {
    //otherwise use consider transition from the continuum
    if(qValueAlpha+energy>0) {
      for(double l=0;l<=maxL_;l+=1.) {
        int piFinal = (int(l)%2==0) ? piInitial : -1*piInitial;
	    for(double jFinal = fabs(l-jInitial); jFinal<=l+jInitial;jFinal+=1.) {
          bool exists = false;
          for(int i = 0;i<spinRatePairs_.size();i++) {
 	     	if(spinRatePairs_[i].Z_==Z-2&&
		       spinRatePairs_[i].A_==A-4&&
		       spinRatePairs_[i].spin_==jFinal&&
		       spinRatePairs_[i].parity_==piFinal) {
		         exists=true;
		         break;
		    }
          }
          if(!exists) {
	    TransitionRateFunc* previous = 
	      (widthCorrectedDecayer_) ? 
	      widthCorrectedDecayer_->spinRatePairs_[spinRatePairs_.size()].rateFunc_ :
	      NULL;
            TransitionRateFunc* newFunc = 
              new TransitionRateFunc(2,4,Z-2,A-4,jInitial,piInitial,
				     jFinal,piFinal,0,1,maxL_,energy,qValueAlpha,
				     totalWidthForCorrection_,uncorrTotalWidthForCorrection_,
				     uncorrTotalWidthSqrdForCorrection_,
				     previous, isCrossSection_);
            spinRatePairs_.push_back(SpinRatePair(Z-2,A-4,jFinal,piFinal,
              qValueAlpha,newFunc,newFunc->Integral()));
            totalIntegral_+=newFunc->Integral();
            totalIntegralSqrd_+=newFunc->Integral()*newFunc->Integral();
	    if(newFunc->GroundStateTransmission()!=0.) 
	      alphaEntrance_ = newFunc->GroundStateTransmission();
	    alphaTotalWidth_+=newFunc->Integral();
          }
        }  
      }
    }
    if(qValueNeutron+energy>0||qValueProton+energy>0) {
      for(double l=0;l<=maxL_;l+=1.) {
	int piFinal = (int(l)%2==0) ? piInitial : -1*piInitial;
	for(double s = fabs(l-0.5); s<=l+0.5; s+=1.) {
	  for(double jFinal = fabs(s-jInitial); jFinal<=s+jInitial;jFinal+=1.) {
	    if(qValueNeutron+energy>0) {
	      bool exists = false;
	      for(int i =0;i<spinRatePairs_.size();i++) {
		if(spinRatePairs_[i].Z_==Z&&
		   spinRatePairs_[i].A_==A-1&&
		   spinRatePairs_[i].spin_==jFinal&&
		   spinRatePairs_[i].parity_==piFinal) {
		  exists=true;
		  break;
		}
	      }
	      if(!exists) {
		TransitionRateFunc* previous = 
		  (widthCorrectedDecayer_) ? 
		  widthCorrectedDecayer_->spinRatePairs_[spinRatePairs_.size()].rateFunc_ :
		  NULL;
		TransitionRateFunc* newFunc = 
		  new TransitionRateFunc(0,1,Z,A-1,jInitial,piInitial,
					 jFinal,piFinal,0.5,1,maxL_,energy,
					 qValueNeutron,totalWidthForCorrection_,uncorrTotalWidthForCorrection_,
					 uncorrTotalWidthSqrdForCorrection_,previous, isCrossSection_);
		spinRatePairs_.push_back(SpinRatePair(Z,A-1,jFinal,piFinal,qValueNeutron,
						      newFunc,newFunc->Integral()));
		totalIntegral_+=newFunc->Integral();
		totalIntegralSqrd_+=newFunc->Integral()*newFunc->Integral();
		if(newFunc->GroundStateTransmission()!=0.) 
		  neutronEntrance_ = newFunc->GroundStateTransmission();
		neutronTotalWidth_+=newFunc->Integral();
	      }
	    }
	    if(qValueProton+energy>0) {
	      bool exists = false;
	      for(int i =0;i<spinRatePairs_.size();i++) {
		if(spinRatePairs_[i].Z_==Z-1&&
		   spinRatePairs_[i].A_==A-1&&
		   spinRatePairs_[i].spin_==jFinal&&
		   spinRatePairs_[i].parity_==piFinal) {
		  exists=true;
		  break;
		}
	      }
	      if(!exists) {
		TransitionRateFunc* previous = 
		  (widthCorrectedDecayer_) ? 
		  widthCorrectedDecayer_->spinRatePairs_[spinRatePairs_.size()].rateFunc_ :
		  NULL;
		TransitionRateFunc* newFunc = 
		  new TransitionRateFunc(1,1,Z-1,A-1,jInitial,piInitial,
					 jFinal,piFinal,0.5,1,maxL_,energy,
					 qValueProton,totalWidthForCorrection_,uncorrTotalWidthForCorrection_,
					 uncorrTotalWidthSqrdForCorrection_,previous,isCrossSection_);
		spinRatePairs_.push_back(SpinRatePair(Z-1,A-1,jFinal,piFinal,qValueProton,
						      newFunc,newFunc->Integral()));
		totalIntegral_+=newFunc->Integral();
		totalIntegralSqrd_+=newFunc->Integral()*newFunc->Integral();
		if(newFunc->GroundStateTransmission()!=0.) 
		  protonEntrance_ = newFunc->GroundStateTransmission();
		protonTotalWidth_+=newFunc->Integral();
	      }
	    }
	  }
	}
      }
    }
    
    for(double jFinal = fabs(jInitial-1.); jFinal<=jInitial+1.;jFinal+=1.) {
      TransitionRateFunc* previous_m1 = 
	(widthCorrectedDecayer_) ? 
	widthCorrectedDecayer_->spinRatePairs_[spinRatePairs_.size()].rateFunc_ :
	NULL;
      TransitionRateFunc* m1Func = 
	new TransitionRateFunc(0,0,Z,A,jInitial,piInitial,jFinal,piInitial,1.,-1,
			       1.,energy,0.,totalWidthForCorrection_,
			       uncorrTotalWidthForCorrection_,uncorrTotalWidthSqrdForCorrection_,
			       previous_m1,isCrossSection_);
      spinRatePairs_.push_back(SpinRatePair(Z,A,jFinal,piInitial,0.,m1Func,m1Func->Integral()));
      totalIntegral_+=m1Func->Integral();
      totalIntegralSqrd_+=m1Func->Integral()*m1Func->Integral();
      gammaEntrance_+= m1Func->GroundStateTransmission();
      gammaTotalWidth_+=m1Func->Integral();
      TransitionRateFunc* previous_e1 = 
	(widthCorrectedDecayer_) ? 
	widthCorrectedDecayer_->spinRatePairs_[spinRatePairs_.size()].rateFunc_ :
	NULL;
      TransitionRateFunc* e1Func = 
	new TransitionRateFunc(0,0,Z,A,jInitial,piInitial,jFinal,-1*piInitial,1.,-1,
			       0.,energy,0.,totalWidthForCorrection_,
			       uncorrTotalWidthForCorrection_,uncorrTotalWidthSqrdForCorrection_,
			       previous_e1,isCrossSection_);
      spinRatePairs_.push_back(SpinRatePair(Z,A,jFinal,-1*piInitial,0.,e1Func,e1Func->Integral()));  	
      totalIntegral_+=e1Func->Integral();
      totalIntegralSqrd_+=e1Func->Integral()*e1Func->Integral();
      gammaEntrance_+= e1Func->GroundStateTransmission();
      gammaTotalWidth_+=e1Func->Integral();
    }
    for(double jFinal = fabs(jInitial-2.); jFinal<=jInitial+2.;jFinal+=1.) {
      TransitionRateFunc* previous_e2 = 
	(widthCorrectedDecayer_) ? 
	widthCorrectedDecayer_->spinRatePairs_[spinRatePairs_.size()].rateFunc_ :
	NULL;
      TransitionRateFunc* e2Func = 
	new TransitionRateFunc(0,0,Z,A,jInitial,piInitial,jFinal,piInitial,1.,-1,
			       2.,energy,0.,totalWidthForCorrection_,
			       uncorrTotalWidthForCorrection_,uncorrTotalWidthSqrdForCorrection_,
			       previous_e2,isCrossSection_);
      spinRatePairs_.push_back(SpinRatePair(Z,A,jFinal,piInitial,0.,e2Func,e2Func->Integral()));
      totalIntegral_+=e2Func->Integral();
      totalIntegralSqrd_+=e2Func->Integral()*e2Func->Integral();
      gammaEntrance_+= e2Func->GroundStateTransmission();
      gammaTotalWidth_+=e2Func->Integral();
    }   
    
    if(!isCrossSection_) BuildCDF();
  }
}

Decayer::~Decayer(){
  for(int i=0 ;i<spinRatePairs_.size();i++) {
    delete spinRatePairs_[i].rateFunc_;
  }
  if(widthCorrectedDecayer_&&widthCorrectedDecayer_->totalWidthForCorrection_!=0.) {
    delete widthCorrectedDecayer_;
  }
}

void Decayer::PrintFunctions() {
  for(std::vector<SpinRatePair>::const_iterator spinRatePair = 
	spinRatePairs_.begin();
      spinRatePair<spinRatePairs_.end();spinRatePair++) {
    if(spinRatePair->rateFunc_==NULL) continue;
    char filename[256];
    char filename2[256];
    char filename3[256];
    char filename4[256];
    sprintf(filename,"levelDensity_Z=%d_A=%d_J=%.1f_Pi=%d.out",spinRatePair->Z_,spinRatePair->A_,
	    spinRatePair->spin_,spinRatePair->parity_);
    sprintf(filename2,"transmissionFunc_Z=%d_A=%d_J=%.1f_Pi=%d.out",spinRatePair->Z_,spinRatePair->A_,
	    spinRatePair->spin_,spinRatePair->parity_);
    sprintf(filename3,"transitionRateFunc_Z=%d_A=%d_J=%.1f_Pi=%d.out",spinRatePair->Z_,spinRatePair->A_,
	    spinRatePair->spin_,spinRatePair->parity_);
    sprintf(filename4,"totalLevelDensity_Z=%d_A=%d.out",spinRatePair->Z_,spinRatePair->A_);
    std::ofstream levelOut(filename);
    std::ofstream transOut(filename2);
    std::ofstream rateOut(filename3);
    std::ofstream totalLevelOut(filename4);
    std::vector<XYPair> function = spinRatePair->rateFunc_->Function();
    for(std::vector<XYPair>::const_iterator it = function.begin();
	it<function.end();it++) {
      levelOut << std::scientific << energy_+spinRatePair->qValue_-it->X_ << ' ' 
	  << spinRatePair->rateFunc_->CalcLevelDensity(energy_+spinRatePair->qValue_-it->X_) 
	  << std::endl;
      totalLevelOut << std::scientific << energy_+spinRatePair->qValue_-it->X_ << ' ' 
	  << spinRatePair->rateFunc_->CalcTotalLevelDensity(energy_+spinRatePair->qValue_-it->X_) 
	  << std::endl;
      transOut << std::scientific << it->X_ << ' ' 
	       << spinRatePair->rateFunc_->CalcTransmissionFunc(it->X_)
	       << std::endl;
      rateOut << std::scientific << it->X_ << ' ' 
	      << it->Y_/totalIntegral_ << ' ' 
	      << std::endl;
    }
    levelOut.flush();
    levelOut.close();
    totalLevelOut.flush();
    totalLevelOut.close();
    transOut.flush();
    transOut.close();
    rateOut.flush();
    rateOut.close();
  }
}

void Decayer::BuildCDF() {
  double offSet=0.;
  for(int pair = 0;pair<spinRatePairs_.size();pair++) {
    std::vector<XYPair> cumulativeSum = 
       spinRatePairs_[pair].rateFunc_->CumulativeSum();
    for(int e = 0;e<cumulativeSum.size();e++) {
      double cdfValue = (cumulativeSum[e].Y_+offSet)/totalIntegral_;
      cdf_.push_back(CDFEntry(pair,cumulativeSum[e].X_,cdfValue));
    }
    offSet+=spinRatePairs_[pair].rateFunc_->Integral();
  }
}

bool Decayer::BuildKnownCDF(int levelIndex, std::vector<Level>& knownLevels) {
  std::vector<GammaTransition> gammas = knownLevels[levelIndex].gammas_;
  for(std::vector<GammaTransition>::const_iterator it = gammas.begin();
      it<gammas.end();it++) {
    totalIntegral_+=it->probability_;
    totalIntegralSqrd_+=it->probability_*it->probability_;
  }
  if(totalIntegral_==0.) return false;
  double offSet =0;
  for(std::vector<GammaTransition>::const_iterator it = gammas.begin();
      it<gammas.end();it++) {
    spinRatePairs_.push_back(SpinRatePair(Z_,A_,knownLevels[it->levelIndex_-1].J_,
					 knownLevels[it->levelIndex_-1].Pi_,0.,
					 NULL,it->probability_));
    double cdfValue = (it->probability_+offSet)/totalIntegral_;
    cdf_.push_back(CDFEntry(spinRatePairs_.size()-1,
			    energy_-knownLevels[it->levelIndex_-1].energy_,
			    cdfValue));
    offSet+=it->probability_;
  }
  return true;
}

void Decayer::PrintCDF() {
  std::ofstream out("cdf.out");
  for(int i=0;i<cdf_.size();i++) {
    int pairIndex = cdf_[i].pairIndex_;
    out << std::fixed << spinRatePairs_[pairIndex].Z_ <<  ' ' << spinRatePairs_[pairIndex].A_ << ' '
	<< spinRatePairs_[pairIndex].spin_ << ' ' << spinRatePairs_[pairIndex].parity_ << ' '
	<< std::scientific << cdf_[i].energy_ << ' ' <<  cdf_[i].value_ << std::endl;
  }
  out.flush();
  out.close();
}

bool Decayer::Decay(int& Z, int& A, double& jFinal, int& piFinal, 
		    double& excitationEnergy, double& decayEnergy) {
  double randomNumber=0.;

  while(randomNumber==0.) randomNumber = double(rand_r(&randomSeed[omp_get_thread_num()]))/double(RAND_MAX);

//MPI_BUILD  while(randomNumber==0.) randomNumber = double(rand())/double(RAND_MAX);

  bool found = false;
  double previousValue = 0.;
  for(std::vector<CDFEntry>::const_iterator it = cdf_.begin();it<cdf_.end();it++) {
    if(previousValue<randomNumber&&randomNumber<=it->value_) {
      Z = spinRatePairs_[it->pairIndex_].Z_;
      A = spinRatePairs_[it->pairIndex_].A_;
      jFinal = spinRatePairs_[it->pairIndex_].spin_;
      piFinal = spinRatePairs_[it->pairIndex_].parity_;
      excitationEnergy = energy_+spinRatePairs_[it->pairIndex_].qValue_-it->energy_;
      decayEnergy = it->energy_;
      found = true;
      break;
    }
  }
  return found;
}

void Decayer::CorrectWidthFluctuations() {
  int maxIt = 2;
  Decayer* previous = this;
  double total = totalIntegral_;
  double uncorrTotal = total; 
  double uncorrTotalSqrd = totalIntegralSqrd_; 
  for(int i = 0;i<maxIt;i++) {
    widthCorrectedDecayer_ = new Decayer(Z_,A_,jInitial_,piInitial_,energy_,
					 total,uncorrTotal,uncorrTotalSqrd,previous);
    previous = widthCorrectedDecayer_;
    total = widthCorrectedDecayer_->totalIntegral_;
  }
}
