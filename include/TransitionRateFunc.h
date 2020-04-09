#ifndef TRANSITIONRATEFUNC_H
#define TRANSITIONRATEFUNC_H

#include <vector>
#include "LevelDensity.h"
#include "TransmissionFunc.h"

class XYPair {
 public:
 XYPair(double X, double Y) : X_(X), Y_(Y) {};
  double X_;
  double Y_;
};

class TransitionRateFunc {
 public:
  TransitionRateFunc(int,int,int,int,
		     double,int,double,int,
		     double,int,double,double,
		     double,double,double,double,
		     TransitionRateFunc*,bool);
  ~TransitionRateFunc() {
    delete transmissionFunc_;
    delete levelDensity_;
  };
  std::vector<XYPair> const Function() {
    return function_;
  };
  std::vector<XYPair> const CumulativeSum() {
    return cumulativeSum_;
  };
  double Integral() const {
    return integral_;
  };
  double CalcLevelDensity(double energy) {
    return levelDensity_->operator()(energy);
  };
  double CalcTransmissionFunc(double energy) {
    return transmissionFunc_->operator()(energy);
  };
  double CalcTotalLevelDensity(double energy) {
    return levelDensity_->TotalLevelDensity(energy);
  };
  double ExclusiveBranching() const {
  	return exclusiveBranching_;
  };
  double GroundStateTransmission() const {
    return groundStateTransmission_;
  };
   static void SetGammaCutoffEnergy(double energy) {
    gammaCutoffEnergy_=energy;
  };
  static double GetGammaCutoffEnergy() {
    return gammaCutoffEnergy_;
  };
 private:
  static const int numCrossSectionSteps_=20;
  static double gammaCutoffEnergy_;
  double integral_;
  double exclusiveBranching_;
  double groundStateTransmission_;
  std::vector<XYPair> function_;
  std::vector<XYPair> cumulativeSum_;
  LevelDensity* levelDensity_;
  TransmissionFunc* transmissionFunc_;
};

#endif