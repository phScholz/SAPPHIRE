#pragma once

#include "TransmissionFunc.h"
#include "NuclearMass.h"
#include "Constants.h"
#include <map>

class SLPair {
 public: 
  SLPair(double s, int l) : s_(s), l_(l) {};
  bool operator<(const SLPair& right) const {
    if(s_ == right.s_ ) return l_ < right.l_;
    else return s_ < right.s_;
  };
  double s_;
  int l_;
};

class ParticleTransmissionFunc : public TransmissionFunc {
 public:
  ParticleTransmissionFunc(int z1, int m1, int z2, int m2, 
			   double jInitial, int piInitial,
			   double jFinal, int piFinal,
			   double spin, int parity, double maxL,
			   double totalWidthForCorrection,
			   double uncorrTotalWidthForCorrection,
			   double uncorrTotalWidthSqrdForCorrection,
			   TransmissionFunc* previous) : 
  TransmissionFunc(z2,m2,jInitial,piInitial,jFinal,piFinal,maxL,
		   totalWidthForCorrection,uncorrTotalWidthForCorrection,
		   uncorrTotalWidthSqrdForCorrection,previous),
    z1_(z1),m1_(m1),parity_(parity),spin_(spin) {
    double mass1, mass2;
    if(!NuclearMass::FindMass(z1_,m1_,mass1)||
       !NuclearMass::FindMass(z2_,m2_,mass2)) redmass_=0.;
    else redmass_=mass1*mass2/(mass1+mass2)/uconv;
    if(z1_==0&&m1_==1) pType_=0;
    else if(z1_==1&&m1_==1) pType_=1;
    else if(z1_==2&&m1_==4) pType_=2;
    else pType_=-1;
  };
  virtual ~ParticleTransmissionFunc() {};
  bool IsValid() {
    if(redmass_==0.||pType_==-1) return false;
    else return true;
  };
  double operator()(double);
  double operator()(double,int);
  void CalcSLDependentFunctions(double,std::map<SLPair,double>&);
  static ParticleTransmissionFunc* 
    CreateParticleTransmissionFunc(int,int,int,int,double,int,
				   double,int,double,int,double,
				   double,double,double,TransmissionFunc*);
  static void SetAlphaFormalism(int formalism) {alphaFormalism_ = formalism;};
  static void SetNeutronFormalism(int formalism) {neutronFormalism_ = formalism;};
  static void SetProtonFormalism(int formalism) {protonFormalism_ = formalism;};
  static void SetPorterThomas(bool);
 protected:
  virtual double CalcTransmission(double,int,double) = 0;
 protected:
  int z1_;
  int m1_;
  int pType_;
  int parity_;
  double redmass_;
  double spin_;
 private:
  static bool porterThomas_;
  static int alphaFormalism_;
  static int protonFormalism_;
  static int neutronFormalism_;
};


