#ifndef PREEQTRANSITIONRATEFUNC_H
#define PREEQTRANSITIONRATEFUNC_H

#include "TransitionRateFunc.h"
#include "ParticleHoleLevelDensity.h"

class PreEqTransitionRateFunc {
 public:
  PreEqTransitionRateFunc(int,int,int,int,int,int,int,
			  int,int,int,int,int,
			  double,int,double,int,
			  double,int,double,double,
			  double);
  ~PreEqTransitionRateFunc() {
  };
  double Integral() const {
    return integral_;
  };
  std::vector<XYPair> const CumulativeSum() {
    return cumulativeSum_;
  };
 private:
  void CalcParticleDecay();
  void CalcNeutronPairProduction();
  void CalcProtonPairProduction();
  void CalcNeutronPairConversion();
  void CalcProtonPairConversion();
 private:
  static const int integrationSteps_ = 100;
  static constexpr double c1_ = 1.;
  static constexpr double c2_ = 1.;
  static constexpr double c3_ = 1.;  
  int z1_;
  int m1_;
  int z2_;
  int m2_;
  int initialNeutronNumber_;
  int initialNeutronHoleNumber_;
  int initialProtonNumber_;
  int initialProtonHoleNumber_;
  int finalNeutronNumber_;
  int finalNeutronHoleNumber_;
  int finalProtonNumber_;
  int finalProtonHoleNumber_;  
  int piInitial_;
  int piFinal_;
  int parity_;
  double jInitial_;
  double jFinal_;
  double spin_;
  double maxL_;
  double compoundE_;
  double qValue_;
  double integral_;
  double M2_;
  std::vector<XYPair> function_;
  std::vector<XYPair> cumulativeSum_;
};

#endif
