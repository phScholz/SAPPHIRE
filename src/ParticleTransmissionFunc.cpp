#include "ParticleTransmissionFunc.h"
#include "EquivSquareWell.h"
#include "McFaddenSatchlerPotential.h"
#include "JLMPotential.h"
#include "Constants.h"
#include <iostream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <omp.h>

extern unsigned int randomSeed[32];


ParticleTransmissionFunc* ParticleTransmissionFunc::CreateParticleTransmissionFunc(int z1, int m1, int z2, int m2, 
							 double jInitial, int piInitial,
							 double jFinal, int piFinal,
							 double spin, int parity, double maxL,
							 double totalWidthForCorrection,
							 double uncorrTotalWidthForCorrection,
							 double uncorrTotalWidthSqrdForCorrection,
							 TransmissionFunc* previous) {
  ParticleTransmissionFunc* transmissionFunc;

  //If Alpha
  if(z1==2&&m1==4) { 	
    if(alphaFormalism_==0) {
      transmissionFunc = new EquivSquareWell(z1,m1,z2,m2, 
					                                  jInitial, piInitial,jFinal,piFinal,
					                                  spin,parity,maxL,
					                                  totalWidthForCorrection,
					                                  uncorrTotalWidthForCorrection,
					                                  uncorrTotalWidthSqrdForCorrection,
					                                  previous);
    } else if(alphaFormalism_==1) {
      transmissionFunc = new McFaddenSatchlerPotential(z1,m1,z2,m2, 
						                                          jInitial, piInitial,jFinal,piFinal,
						                                          spin,parity,maxL,
						                                          totalWidthForCorrection,
						                                          uncorrTotalWidthForCorrection,
						                                          uncorrTotalWidthSqrdForCorrection,
						                                          previous);
    } else {
      std::cout << "Specified alpha transmission formalism doesn't exist.  Exiting." << std::endl;
      exit(1);	
    }
  }else if((z1==1&&m1==1&&protonFormalism_==1) || (z1==0&&m1==1&&neutronFormalism_==1)) {
      transmissionFunc = new JLMPotential(z1,m1,z2,m2, 
					                                jInitial, piInitial,jFinal,piFinal,
					                                spin,parity,maxL,
					                                totalWidthForCorrection,
					                                uncorrTotalWidthForCorrection,
					                                uncorrTotalWidthSqrdForCorrection,
					                                previous);
  } else if ((z1==1&&m1==1&&protonFormalism_==0) || (z1==0&&m1==1&&neutronFormalism_==0)) {
      transmissionFunc = new EquivSquareWell(z1,m1,z2,m2, 
					                                jInitial, piInitial,jFinal,piFinal,
					                                spin,parity,maxL,
					                                totalWidthForCorrection,
					                                uncorrTotalWidthForCorrection,
					                                uncorrTotalWidthSqrdForCorrection,
					                                previous);
  } else {
      std::cout << "Specified proton/neutron transmission formalism doesn't exist.  Exiting." << std::endl;
      exit(1);	
  }
  
  return transmissionFunc;
}


double ParticleTransmissionFunc::operator()(double energy){
  gsl_rng *r;
  if(porterThomas_) {
    const gsl_rng_type *T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    gsl_rng_set (r, rand_r(&randomSeed[omp_get_thread_num()])*time(NULL));
  }

  std::map<SLPair,double> functions;
  CalcSLDependentFunctions(energy,functions);
  double sum=0.;
  int i=0;

  for(std::map<SLPair,double>::const_iterator it = functions.begin(); it!=functions.end();it++) {
      double chirand = (porterThomas_) ? gsl_ran_chisq (r, 1.) : 1.;
      
      if(totalWidthForCorrection_!=0.) {
        double previousT = (previous_) ? ((ParticleTransmissionFunc*)previous_)->operator()(energy,i) : it->second;
        double tBar = uncorrTotalWidthSqrdForCorrection_/uncorrTotalWidthForCorrection_;
        double exponent = 4.*tBar/uncorrTotalWidthForCorrection_*(1.+it->second/uncorrTotalWidthForCorrection_)/(1.+3*tBar/uncorrTotalWidthForCorrection_);
        double WCF = 1.+2./(1.+pow(it->second,exponent))+87.*pow((it->second-tBar)/uncorrTotalWidthForCorrection_,2.)*pow(it->second/uncorrTotalWidthForCorrection_,5.);
        
        sum+=it->second/(1.+(WCF-1.)*previousT/totalWidthForCorrection_);
    
    } else sum+=it->second*chirand;
      i++;
  }
  if(porterThomas_) gsl_rng_free (r);
  return sum;
}

double ParticleTransmissionFunc::operator()(double energy, int which) {
  std::map<SLPair,double> functions;
  CalcSLDependentFunctions(energy,functions);
  double sum=0.;
  int i=0;
  for(std::map<SLPair,double>::const_iterator it = functions.begin();
      it!=functions.end();it++) {
    if(i==which) {
      if(totalWidthForCorrection_!=0.) {
	double previousT = (previous_) ? previous_->operator()(energy) : it->second;
	double tBar = uncorrTotalWidthSqrdForCorrection_/uncorrTotalWidthForCorrection_;
	double exponent = 4.*tBar/uncorrTotalWidthForCorrection_
	  *(1.+it->second/uncorrTotalWidthForCorrection_)/(1.+3*tBar/uncorrTotalWidthForCorrection_);
	double WCF = 1.+2./(1.+pow(it->second,exponent))+87.*pow((it->second-tBar)/uncorrTotalWidthForCorrection_,2.)*
	  pow(it->second/uncorrTotalWidthForCorrection_,5.);
	sum=it->second/(1.+(WCF-1.)*previousT/totalWidthForCorrection_);
      } else sum=it->second;
      break;
    }
    i++;
  }
  return sum;
}

void ParticleTransmissionFunc::CalcSLDependentFunctions(double energy,
							std::map<SLPair,double>& functions) {
  for(double s=fabs(jFinal_-spin_);s<=jFinal_+spin_;s+=1.0) {
    for(double l=fabs(jInitial_-s);l<=jInitial_+s;l+=1.0) {
      double intPart;
      if(modf(l,&intPart)!=0.) {
	std::cout << "WARNING: Calculated L not an integer." << std::endl;
	continue;
      }
      if(int(intPart)>maxL_) break;
      int parity = (int(intPart)%2==0) ? parity_*piFinal_ : -1*parity_*piFinal_;
      if(piInitial_!=parity) continue;
      functions[SLPair(s,int(intPart))]=CalcTransmission(s,int(intPart),energy);
    }
  }
}
