#include "Decayer/Decayer.h"
#include "TransitionRateFunc.h"
#include "Databases/NuclearMass.h"
#include "Databases/NuclearLevels.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

extern unsigned int randomSeed[32];

bool Decayer::BoundStateCheck(){
  bool knownCDFBuilt = false;

  //Looking for known levels of the nucleus
  std::vector<Level> knownLevels = NuclearLevels::FindLevels(Z_,A_);

  //Check if the state is a known state ... accuracy +/- 1.e-10 MeV and store the level in foundLevel
  //Otherwise foundLevel will be the last known level
  std::vector<Level>::const_iterator foundLevel = knownLevels.end();
  if(qValueAlpha_+energy_<=0.&& qValueProton_+energy_<=0.&& qValueNeutron_+energy_<=0.) {
    for(std::vector<Level>::const_iterator it = knownLevels.begin(); it<knownLevels.end();++it) {
      if((it->energy_-1.e-10<=energy_ && energy_ <=it->energy_+1.e-10) && jInitial_==it->J_&&piInitial_==it->Pi_){
	      foundLevel = it;
	      break;
      }
    }  
  }

  if(foundLevel!=knownLevels.end() && foundLevel->gammas_.size()>0) {
    //if level is known, and has a known decay scheme use it
    knownCDFBuilt = BuildKnownCDF(int(foundLevel-knownLevels.begin()),knownLevels);
  } 

  return knownCDFBuilt;
}

void Decayer::InitializeQValues(){
    if(!NuclearMass::QValue(Z_,A_,Z_,A_-1,qValueNeutron_)) {
    std::cout << "Unknown masses requested.  Aborting." << std::endl;
    exit(1);
  }
  if(!NuclearMass::QValue(Z_,A_,Z_-1,A_-1,qValueProton_)) {
    std::cout << "Unknown masses requested.  Aborting." << std::endl;
    exit(1);
  }
  if(!NuclearMass::QValue(Z_,A_,Z_-2,A_-4,qValueAlpha_)) {
    std::cout << "Unknown masses requested.  Aborting." << std::endl;
    exit(1);
  }
}

void Decayer::InitializeWidths(){
  neutronEntrance_=0.;
  gammaEntrance_=0.;
  protonEntrance_=0.;
  alphaEntrance_=0.;
  neutronTotalWidth_=0.;
  alphaTotalWidth_=0.;
  gammaTotalWidth_=0.;
  protonTotalWidth_=0.;

  //total integrals and the square are set to 0
  totalIntegral_=0.;    
  totalIntegralSqrd_=0.; 
}

void Decayer::SetAlphaSpinRatePairs(){
  if(alphaDecay_){
    for(double l=0;l<=maxL_;l+=1.) {
        int piFinal = (int(l)%2==0) ? piInitial_ : -1*piInitial_;

	      for(double jFinal = fabs(l-jInitial_); jFinal<=l+jInitial_;jFinal+=1.) {
          bool exists = false;
          
          //Check if this spin rate pair already exists or not
          for(int i = 0;i<spinRatePairs_.size();i++) {
 	     	    if(spinRatePairs_[i].Z_==Z_-2&&
		          spinRatePairs_[i].A_==A_-4&&
		          spinRatePairs_[i].spin_==jFinal&&
		          spinRatePairs_[i].parity_==piFinal) {
		          exists=true;
		          break;
		        }
          }
          
          //If this pair doesn't exist yet, create it
          if(!exists) {
	            TransitionRateFunc* previous = (widthCorrectedDecayer_) ? widthCorrectedDecayer_->spinRatePairs_[spinRatePairs_.size()].rateFunc_ : NULL;
              TransitionRateFunc* newFunc = new TransitionRateFunc(2,4,Z_-2,A_-4,jInitial_,piInitial_, jFinal,piFinal,0,1,maxL_,energy_,qValueAlpha_, totalWidthForCorrection_,uncorrTotalWidthForCorrection_, uncorrTotalWidthSqrdForCorrection_, previous, isCrossSection_);
              spinRatePairs_.push_back(SpinRatePair(Z_-2,A_-4,jFinal,piFinal, qValueAlpha_,newFunc,newFunc->Integral()));
              totalIntegral_+=newFunc->Integral();
              totalIntegralSqrd_+=newFunc->Integral()*newFunc->Integral();
	    
              if(newFunc->GroundStateTransmission()!=0.) 
	              alphaEntrance_ = newFunc->GroundStateTransmission();
	              alphaTotalWidth_+=newFunc->Integral();
          }
        }  
    }
  }
}

void Decayer::SetProtonNeutronSpinRatePairs(){
        for(double l=0;l<=maxL_;l+=1.) {
	          int piFinal = (int(l)%2==0) ? piInitial_ : -1*piInitial_;
	          for(double s = fabs(l-0.5); s<=l+0.5; s+=1.) {
	            for(double jFinal = fabs(s-jInitial_); jFinal<=s+jInitial_;jFinal+=1.) {
	              //if neutron decay is  energetically possible
                if(qValueNeutron_+energy_>0 && neutronDecay_) {
	                bool exists = false;
	                  for(int i =0;i<spinRatePairs_.size();i++) {
		                  if(spinRatePairs_[i].Z_==Z_&&
		                     spinRatePairs_[i].A_==A_-1&&
		                     spinRatePairs_[i].spin_==jFinal&&
		                     spinRatePairs_[i].parity_==piFinal) {
		                    exists=true;
		                    break;
		                  }
	                  }
	              
                  if(!exists) {
		                TransitionRateFunc* previous = (widthCorrectedDecayer_) ? widthCorrectedDecayer_->spinRatePairs_[spinRatePairs_.size()].rateFunc_ : NULL;
		                TransitionRateFunc* newFunc =  new TransitionRateFunc(0,1,Z_,A_-1,jInitial_,piInitial_, jFinal,piFinal,0.5,1,maxL_,energy_, qValueNeutron_,totalWidthForCorrection_,uncorrTotalWidthForCorrection_, uncorrTotalWidthSqrdForCorrection_,previous, isCrossSection_);
                    spinRatePairs_.push_back(SpinRatePair(Z_,A_-1,jFinal,piFinal,qValueNeutron_, newFunc,newFunc->Integral()));
		                totalIntegral_+=newFunc->Integral();
		                totalIntegralSqrd_+=newFunc->Integral()*newFunc->Integral();
		              
                    if(newFunc->GroundStateTransmission()!=0.) neutronEntrance_ = newFunc->GroundStateTransmission();
		              
                    neutronTotalWidth_+=newFunc->Integral();
	                }
	              }
                //If proton decay is energetically possible
                if(qValueProton_+energy_>0 && protonDecay_) {
	                bool exists = false;
	                for(int i =0;i<spinRatePairs_.size();i++) {
		                if(spinRatePairs_[i].Z_==Z_-1&&
		                  spinRatePairs_[i].A_==A_-1&&
		                  spinRatePairs_[i].spin_==jFinal&&
		                  spinRatePairs_[i].parity_==piFinal) {
		                  exists=true;
		                  break;
		                }
	                }
	                
                  if(!exists) {
		                TransitionRateFunc* previous = (widthCorrectedDecayer_) ? widthCorrectedDecayer_->spinRatePairs_[spinRatePairs_.size()].rateFunc_ : NULL;
		                TransitionRateFunc* newFunc = new TransitionRateFunc(1,1,Z_-1,A_-1,jInitial_,piInitial_, jFinal,piFinal,0.5,1,maxL_,energy_, qValueProton_,totalWidthForCorrection_,uncorrTotalWidthForCorrection_, uncorrTotalWidthSqrdForCorrection_,previous,isCrossSection_);
		                spinRatePairs_.push_back(SpinRatePair(Z_-1,A_-1,jFinal,piFinal,qValueProton_,
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

void Decayer::SetE1M1SpinRatePairs(){
  if(gammaDecay_){
    for(double jFinal = fabs(jInitial_-1.); jFinal<=jInitial_+1.;jFinal+=1.) {
      //Define transitionRateFuncs for M1
      TransitionRateFunc* previous_m1 = (widthCorrectedDecayer_) ? widthCorrectedDecayer_->spinRatePairs_[spinRatePairs_.size()].rateFunc_ : NULL;
      TransitionRateFunc* m1Func = 	new TransitionRateFunc(0,0,Z_,A_,jInitial_,piInitial_,jFinal,piInitial_,1.,-1, 1.,energy_,0.,totalWidthForCorrection_, uncorrTotalWidthForCorrection_,uncorrTotalWidthSqrdForCorrection_, previous_m1,isCrossSection_);
      
      //Push back allowed transitions
      spinRatePairs_.push_back(SpinRatePair(Z_,A_,jFinal,piInitial_,0.,m1Func,m1Func->Integral()));
      totalIntegral_+=m1Func->Integral();
      totalIntegralSqrd_+= m1Func->Integral()*m1Func->Integral();
      gammaEntrance_+= m1Func->GroundStateTransmission();
      gammaTotalWidth_+=m1Func->Integral();

      //Define transitionRateFuncs for E1
      TransitionRateFunc* previous_e1 = (widthCorrectedDecayer_) ? widthCorrectedDecayer_->spinRatePairs_[spinRatePairs_.size()].rateFunc_ : NULL;
      TransitionRateFunc* e1Func = new TransitionRateFunc(0,0,Z_,A_,jInitial_,piInitial_,jFinal,-1*piInitial_,1.,-1, 0.,energy_,0.,totalWidthForCorrection_, uncorrTotalWidthForCorrection_,uncorrTotalWidthSqrdForCorrection_, previous_e1,isCrossSection_);
      
      //Push back allowed transitions
      spinRatePairs_.push_back(SpinRatePair(Z_,A_,jFinal,-1*piInitial_,0.,e1Func,e1Func->Integral()));  	
      totalIntegral_+=e1Func->Integral();
      totalIntegralSqrd_+=e1Func->Integral()*e1Func->Integral();
      gammaEntrance_+= e1Func->GroundStateTransmission();
      gammaTotalWidth_+=e1Func->Integral();
    }
  }
}

void Decayer::SetE2SpinRatePairs(){
  if(gammaDecay_){
    for(double jFinal = fabs(jInitial_-2.); jFinal<=jInitial_+2.;jFinal+=1.) {
      //Define TransitionRateFuncs for E2
      TransitionRateFunc* previous_e2 =	(widthCorrectedDecayer_) ? widthCorrectedDecayer_->spinRatePairs_[spinRatePairs_.size()].rateFunc_ : NULL;
      TransitionRateFunc* e2Func = new TransitionRateFunc(0,0,Z_,A_,jInitial_,piInitial_,jFinal,piInitial_,1.,-1, 2.,energy_,0.,totalWidthForCorrection_, uncorrTotalWidthForCorrection_,uncorrTotalWidthSqrdForCorrection_, previous_e2,isCrossSection_);
      
      //Push back allowed transitions
      spinRatePairs_.push_back(SpinRatePair(Z_,A_,jFinal,piInitial_,0.,e2Func,e2Func->Integral()));
      totalIntegral_+=e2Func->Integral();
      totalIntegralSqrd_+=e2Func->Integral()*e2Func->Integral();
      gammaEntrance_+= e2Func->GroundStateTransmission();
      gammaTotalWidth_+=e2Func->Integral();
    }
  }
}

Decayer::Decayer(int Z, int A, double jInitial, int piInitial, double energy, double totalWidthForCorrection, double uncorrTotalWidthForCorrection, double uncorrTotalWidthSqrdForCorrection,		 Decayer* widthCorrectedDecayer) :
  Z_(Z), A_(A), piInitial_(piInitial), jInitial_(jInitial), energy_(energy), totalWidthForCorrection_(totalWidthForCorrection), uncorrTotalWidthForCorrection_(uncorrTotalWidthForCorrection), uncorrTotalWidthSqrdForCorrection_(uncorrTotalWidthSqrdForCorrection), widthCorrectedDecayer_(widthCorrectedDecayer) {

  InitializeWidths();

  InitializeQValues();

  //Check whether it is a bound state, then construct a knownCDF
  if(!BoundStateCheck()) {
    //otherwise  consider transition from the continuum
    
    //qValueAlpha_=0;
    //If Alphadecay is energetically possible
    if(qValueAlpha_+energy_>0) {
      SetAlphaSpinRatePairs();
    }

    //qValueNeutron_=0;
    //qValueProton_=0;
    //If neutron OR proton decay is energetically possible
    if(qValueNeutron_+energy_>0||qValueProton_+energy>0) {
      SetProtonNeutronSpinRatePairs();    
    }  
    
    //This is now for gammas with L=1 -> E1, M1
    SetE1M1SpinRatePairs();

    //This is now for gammas with L=2 -> E2, M2
    SetE2SpinRatePairs();
    
    //Build the cummulative distribution function, if this is not a cross section calculation
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
  for(std::vector<SpinRatePair>::const_iterator spinRatePair = spinRatePairs_.begin(); spinRatePair<spinRatePairs_.end(); ++spinRatePair) {
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
    
    for(std::vector<XYPair>::const_iterator it = function.begin(); it<function.end(); ++it) {
      levelOut << std::scientific << energy_+spinRatePair->qValue_-it->X_ << ' ' << spinRatePair->rateFunc_->CalcLevelDensity(energy_+spinRatePair->qValue_-it->X_)  << std::endl;
      totalLevelOut << std::scientific << energy_+spinRatePair->qValue_-it->X_ << ' '  << spinRatePair->rateFunc_->CalcTotalLevelDensity(energy_+spinRatePair->qValue_-it->X_) << std::endl;
      transOut << std::scientific << it->X_ << ' ' << spinRatePair->rateFunc_->CalcTransmissionFunc(it->X_) << std::endl;
      rateOut << std::scientific << it->X_ << ' '  << it->Y_/totalIntegral_ << ' ' << std::endl;
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
  for(int pair = 0; pair<spinRatePairs_.size(); pair++) {
    std::vector<XYPair> cumulativeSum = spinRatePairs_[pair].rateFunc_->CumulativeSum();
    for(int e = 0;e<cumulativeSum.size();e++) {
      double cdfValue = (cumulativeSum[e].Y_+offSet)/totalIntegral_;
      cdf_.push_back(CDFEntry(pair,cumulativeSum[e].X_,cdfValue));
    }
    offSet+=spinRatePairs_[pair].rateFunc_->Integral();
  }
}

bool Decayer::BuildKnownCDF(int levelIndex, std::vector<Level>& knownLevels) {
  //Create a vector for all the gamma transitions known for the specific level.
  std::vector<GammaTransition> gammas = knownLevels[levelIndex].gammas_;

  //For loop to calculate the totaldecay width and its squared for all the known gamma transitions of the level
  for(std::vector<GammaTransition>::const_iterator it = gammas.begin(); it<gammas.end(); ++it) {
    totalIntegral_+=it->probability_;
    totalIntegralSqrd_+=it->probability_*it->probability_;
  }

  //Return false if the totalwidth is 0.
  if(totalIntegral_==0.) return false;

  //Some kind of offSet for the partial width
  double offSet =0;

  //For all the known gammatransitions, a spinRatePair is created and pushed back into the spinRatePairs_
  //The cdfValue is calculated as the ratio of the partial decay probability and the total widht
  //A CDFEnty is generated for the last spinRatePair using the number of the last entry of SpinRatePairs, the energy difference and the cdfValue
  //Since its a cummulative distribution function the probability of the previous entry is the offset of the next
  for(std::vector<GammaTransition>::const_iterator it = gammas.begin(); it<gammas.end(); ++it) {
    spinRatePairs_.push_back(SpinRatePair(Z_,A_,knownLevels[it->levelIndex_-1].J_, knownLevels[it->levelIndex_-1].Pi_,0., NULL, it->probability_));
    double cdfValue = (it->probability_+offSet)/totalIntegral_;
    cdf_.push_back(CDFEntry(spinRatePairs_.size()-1, energy_-knownLevels[it->levelIndex_-1].energy_, cdfValue));
    offSet+=it->probability_;
  }
  return true;
}

void Decayer::PrintCDF() {
  std::ofstream out("cdf.out");
  for(int i=0;i<cdf_.size();i++) {
    int pairIndex = cdf_[i].pairIndex_;
    out << std::fixed << spinRatePairs_[pairIndex].Z_ <<  ' ' << spinRatePairs_[pairIndex].A_ << ' '	<< spinRatePairs_[pairIndex].spin_ << ' ' << spinRatePairs_[pairIndex].parity_ << ' '	<< std::scientific << cdf_[i].energy_ << ' ' <<  cdf_[i].value_ << std::endl;
  }
  out.flush();
  out.close();
}

bool Decayer::Decay(int& Z, int& A, double& jFinal, int& piFinal, 
		    double& excitationEnergy, double& decayEnergy) {
  double randomNumber=0.;
  while(randomNumber==0.) randomNumber = double(rand_r(&randomSeed[omp_get_thread_num()]))/double(RAND_MAX);

  bool found = false;
  
  for(std::vector<CDFEntry>::const_iterator it = cdf_.begin(); it<cdf_.end(); ++it) {
    if(randomNumber<=it->value_) {
      Z = spinRatePairs_[it->pairIndex_].Z_;
      A = spinRatePairs_[it->pairIndex_].A_;
      jFinal = spinRatePairs_[it->pairIndex_].spin_;
      piFinal = spinRatePairs_[it->pairIndex_].parity_;
      excitationEnergy = energy_ + spinRatePairs_[it->pairIndex_].qValue_ - it->energy_;
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
    widthCorrectedDecayer_ = new Decayer(Z_,A_,jInitial_,piInitial_,energy_, total,uncorrTotal,uncorrTotalSqrd,previous);
    previous = widthCorrectedDecayer_;    
    total = widthCorrectedDecayer_->totalIntegral_;
  }
}
