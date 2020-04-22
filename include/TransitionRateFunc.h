#pragma once

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
 /**
  * @brief Constructor of a TransitionRateFunc
  * @param z1
  * @param m1
  * @param z2
  * @param m2 
  * @param jInitial
  * @param piInitial  
  * @param jFinal 
  * @param piFinal
  * @param spin
  * @param parity
  * @param maxL
  * @param compoundE
  * @param qValue
  * @param totalWidthForCorrection
  * @param uncorrTotalWidthForCorrection
  * @param uncorrTotalWidthSqrdForCorrection
  * @param previous
  * @param isCrossSection
  * @details
  * 1. Choose level density model (only Rauscher atm)
  * 2. Initialize a transmissionFunc_ depending on the projectile
  * 3. Read in levels of nucleus via NuclearLevels::FindLevels()
  * 4. 
  */
  TransitionRateFunc(int z1, int m1, int z2, int m2,
				             double jInitial, int piInitial, 
				             double jFinal, int piFinal,
				             double spin, int parity, double maxL,
				             double compoundE, double qValue,
				             double totalWidthForCorrection, 
				             double uncorrTotalWidthForCorrection, 
				             double uncorrTotalWidthSqrdForCorrection, 
				             TransitionRateFunc* previous,
				             bool isCrossSection);
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
  LevelDensity* levelDensity_;          /**< LevelDensity model*/
  TransmissionFunc* transmissionFunc_; /**< Particle or GammaTransmissionFunc*/
};

