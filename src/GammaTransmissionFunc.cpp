#include "GammaTransmissionFunc.h"
#include "BrinkAxelGSF.h"
#include "KopeckyUhlGSF.h"
#include "McCullaghGSF.h"
#include "Constants.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <omp.h>



extern unsigned int randomSeed[12];


void GammaTransmissionFunc::SetPorterThomas(bool yn) {
  porterThomas_=yn;
}

void GammaTransmissionFunc::SetEGDRType(int type) {
  egdrType_=type;
}

void GammaTransmissionFunc::InitializeGDRParameters(std::string filename) {
  std::cout << "Reading GDR parameters file..." << std::endl;
  std::ifstream in(filename.c_str());
  if(!in) {
    std::cout << "Could not read GDR parameters file." << std::endl;
    exit(1);
  }
  std::string line;
  for(int i = 0;i<4;++i) std::getline(in,line);
  while(!in.eof()) {
    std::getline(in,line);
    if(in.eof()) continue;
    std::istringstream lineStream(line);
    int Z,A; lineStream>>Z>>A;
    lineStream.ignore(3);
    double eta,e1,w1,e2,w2;
    lineStream>>eta>>e1>>w1>>e2>>w2;
    gdrTable_[MassKey(Z,A)]=GDRParameters(eta,e1,w1,0.,e2,w2,0.);
  }
  in.close();
}

GammaTransmissionFunc::GammaTransmissionFunc(int z2, int m2, double jInitial, int piInitial,
					     double jFinal, int piFinal, double maxL, 
					     double totalWidthForCorrection,
					     double uncorrTotalWidthForCorrection,
					     double uncorrTotalWidthSqrdForCorrection,
					     TransmissionFunc* previous) : 
  TransmissionFunc(z2,m2,jInitial,piInitial,jFinal,piFinal,maxL,
		   totalWidthForCorrection,uncorrTotalWidthForCorrection,
		   uncorrTotalWidthSqrdForCorrection,previous) {
  double k =1./pi/pi/hbarc/hbarc/10.;
  if(maxL==0.) {
    k/=3.;
    double protonMass;
    NuclearMass::FindMass(1,1,protonMass);
    double kSigmaGamma = k*4.*pi*hbarc*hbarc*10.*fstruc*(m2_-z2_)*z2_/m2_*1.2/protonMass;
    GDRTable::const_iterator parameters = gdrTable_.find(MassKey(z2,m2));
    if(parameters!=gdrTable_.end()) {
      GDRParameters tempParameters = parameters->second;
      tempParameters.kSigmaGamma_[0] = 1./3.*kSigmaGamma;
      tempParameters.kSigmaGamma_[1] = 2./3.*kSigmaGamma;
      gdrParameters_ = tempParameters;
    } else {
      gdrParameters_ = GDRParameters(1.,
				     31.2*pow(double(m2_),-0.33333333333)+
				     20.6*pow(double(m2_),-0.16666666667),
				     0.026*pow(31.2*pow(double(m2_),-0.33333333333)+
					       20.6*pow(double(m2_),-0.16666666667),1.91),
				     kSigmaGamma,
				     1.,0.,0.);
    }
  } else if(maxL==1.) {
    k/=3.;
    double kSigmaGamma = 1.58e-9*pow(double(m2_),0.47)*
      (pow(49.-pow(41.*pow(double(m2_),-0.33333333333),2.),2.)+784.)/28.;
    double GDREnergy = 41.*pow(double(m2_),-0.33333333333);
    double GDRWidth = 4.;
    gdrParameters_ = GDRParameters(1.,
				   GDREnergy,
				   GDRWidth,
				   kSigmaGamma,
				   1.,0.,0.);
  } else if(maxL==2.) {
    k/=5.;
    double kSigmaGamma = k*0.59535*z2_*z2_/double(m2_);
    double GDREnergy = 63.*pow(double(m2_),-0.33333333333);
    double GDRWidth = 6.11-0.012*m2_;
    gdrParameters_ = GDRParameters(1.,
				   GDREnergy,
				   GDRWidth,
				   kSigmaGamma,
				   1,0.,0.);
  } 
}

GammaTransmissionFunc*
GammaTransmissionFunc::CreateGammaTransmissionFunc(int z2, int m2, double jInitial, int piInitial,
						   double jFinal, int piFinal, double maxL, 
						   LevelDensityFormula* levelDensity, double totalWidthForCorrection,
						   double uncorrTotalWidthForCorrection,
						   double uncorrTotalWidthSqrdForCorrection,
						   TransmissionFunc* previous, double compoundE) {
  //Variable maxL determines multipolarity (E1 = 0., M1= 1., E2 = 2.).

  GammaTransmissionFunc* function;
  if((maxL==0.&&egdrType_==0)||
     (maxL==1.&&mgdrType_==0)|| 
     (maxL==2.&&egqrType_==0)) function = new BrinkAxelGSF(z2,m2,jInitial,piInitial,jFinal,piFinal,maxL, 
							  totalWidthForCorrection,
							  uncorrTotalWidthForCorrection,
							  uncorrTotalWidthSqrdForCorrection,
							  previous);
  else if(maxL==0.&&egdrType_==1) function = new KopeckyUhlGSF(z2,m2,jInitial,piInitial,jFinal,piFinal,maxL, 
							      totalWidthForCorrection,
							      uncorrTotalWidthForCorrection,
							      uncorrTotalWidthSqrdForCorrection,
							      previous,levelDensity,compoundE);
  else if((maxL==0.&&egdrType_==2)||
	  (maxL==1.&&mgdrType_==2)|| 
	  (maxL==2.&&egqrType_==2)) function = new McCullaghGSF(z2,m2,jInitial,piInitial,jFinal,piFinal,maxL, 
							       totalWidthForCorrection,
							       uncorrTotalWidthForCorrection,
							       uncorrTotalWidthSqrdForCorrection,
							       previous);
  else {
    std::cout << "Giant resonance shape / multipolarity combination not known.  Aborting." << std::endl;
    std::exit(1);
  }
  return function;
}
 
double GammaTransmissionFunc::operator()(double energy) {
  double chirand = 1.;
  if(porterThomas_) {
    const gsl_rng_type *Tr = gsl_rng_default;
    gsl_rng *r = gsl_rng_alloc (Tr);

    gsl_rng_set (r, rand_r(&randomSeed[omp_get_thread_num()])*time(NULL));
    chirand = gsl_ran_chisq (r, 1.);
    gsl_rng_free (r);
  }
  double T = (maxL_==0.||maxL_==1.) ? 2.*pi*CalcStrengthFunction(energy)*pow(energy,3.) : 
    2.*pi*CalcStrengthFunction(energy)*pow(energy,5.);
  if(totalWidthForCorrection_!=0.) {
    double previousT = (previous_) ? previous_->operator()(energy) : T;
    double tBar = uncorrTotalWidthSqrdForCorrection_/uncorrTotalWidthForCorrection_;
    double exponent = 4.*tBar/uncorrTotalWidthForCorrection_
      *(1.+T/uncorrTotalWidthForCorrection_)/(1.+3*tBar/uncorrTotalWidthForCorrection_);
    double WCF = 1.+2./(1.+pow(T,exponent))+87.*pow((T-tBar)/uncorrTotalWidthForCorrection_,2.)*
      pow(T/uncorrTotalWidthForCorrection_,5.);
    return T/(1.+(WCF-1.)*previousT/totalWidthForCorrection_);
  } else return T*chirand;
}


