/**
 * @file TransitionRateFunc.h
 * @brief Declaration of the base class for the calculation of transition rates.
 * @date 2020-04-27
 */
#pragma once

#include <vector>
#include "LevelDensityFormula.h"
#include "TransmissionFunc.h"

/**
 * @brief Simple pair of two double variables (x,y)
 */
class XYPair {
 public:
 XYPair(double X, double Y) : X_(X), Y_(Y) {};
  double X_;
  double Y_;
};

/**
 * @brief Class to calculate the transition rates to different energy bins or levels
 */
class TransitionRateFunc {
 public:
 /**
  * @brief Constructor of a TransitionRateFunc
  * @param z1 charge of ejectile
  * @param m1 mass of ejectile
  * @param z2 charge of resiudal
  * @param m2 mass of residual
  * @param jInitial spin of compound
  * @param piInitial  parity of compound
  * @param jFinal spin of residual
  * @param piFinal parity of residual
  * @param spin spin of ejectile
  * @param parity parity of ejectile
  * @param maxL maximum angular momentum
  * @param compoundE Excitation energy before decay.
  * @param qValue QValue for the decay.
  * @param totalWidthForCorrection
  * @param uncorrTotalWidthForCorrection
  * @param uncorrTotalWidthSqrdForCorrection
  * @param previous
  * @param isCrossSection
  * @details
  * 1. Choose level density model (only Rauscher atm)
  * 2. Initialize a transmissionFunc_ depending on the projectile
  * 3. Read in levels of nucleus via NuclearLevels::FindLevels()
  * 4. Now rates will be calculated and integrated in two steps:
  *   - Contiuum loop: \f[ \int r(E\prime) dE\prime  = \int \rho(E_x-E\prime, J, \Pi) \cdot T(E\prime, J, \Pi)\f]. 
  *     This integral is approximated by the [composite simpson rule](https://en.wikipedia.org/wiki/Simpson%27s_rule#Composite_Simpson's_rule).
  *   - Discrete loop: \f[ \sum_i r_i = \sum_i T(E_i, J, \Pi) \f]
  *     The sum runs over all known levels of the specific jFinal and piFinal.
  * 5. The sum of both integrals is stored in integral_; 
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
  
  /**
   * @brief Deconstructor
   */
  ~TransitionRateFunc() {
    delete transmissionFunc_;
    delete levelDensity_;
  };

  /** Getter Function**/
  std::vector<XYPair> const Function() {
    return function_;
  };

  /** Getter cumulativeSum_**/
  std::vector<XYPair> const CumulativeSum() {
    return cumulativeSum_;
  };

  /** Getter integral**/
  double Integral() const {
    return integral_;
  };

  /** Getter leveldensity at excitation energy**/
  double CalcLevelDensity(double energy) {
    return levelDensity_->operator()(energy);
  };

  /** Getter transmission coefficient for specific energy**/
  double CalcTransmissionFunc(double energy) {
    return transmissionFunc_->operator()(energy);
  };

  /** Getter total level density at energy**/
  double CalcTotalLevelDensity(double energy) {
    return levelDensity_->TotalLevelDensity(energy);
  };

  /** Getter for exclusive branching**/
  double ExclusiveBranching() const {
  	return exclusiveBranching_;
  };

  /** Getter groundStateTransmission**/
  double GroundStateTransmission() const {
    return groundStateTransmission_;
  };

  /** Setter gammaCutOff**/
   static void SetGammaCutoffEnergy(double energy) {
    gammaCutoffEnergy_=energy;
  };

  /** Getter gammaCutOff**/
  static double GetGammaCutoffEnergy() {
    return gammaCutoffEnergy_;
  };
  
 private:
  static const int numCrossSectionSteps_=20;
  static double gammaCutoffEnergy_; /**< Excitation energy until which the gamma-decay is assumed to be exclusive.*/
  double integral_;
  double exclusiveBranching_;
  double groundStateTransmission_;
  std::vector<XYPair> function_;
  std::vector<XYPair> cumulativeSum_;
  LevelDensityFormula* levelDensity_;          /**< LevelDensity model*/
  TransmissionFunc* transmissionFunc_; /**< Particle or GammaTransmissionFunc*/
};

