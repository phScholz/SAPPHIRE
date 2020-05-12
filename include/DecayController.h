/**
 * @file DecayController.h
 * @brief Declaration of class DecayController
 */

#pragma once

#include <vector>
#include "Decayer.h"
#include "DecayProduct.h"

class TVector3;

/**
 * @brief Class to controll a decay.
 */
class DecayController {
 
 public:
  /**
   * @brief Constructor. Initializes members.
   * @param Z Charge number
   * @param A Mass number
   * @param jInitial Spin
   * @param piInitial, parity
   * @param energy Excitation energy
   * @param initialNeutronNumber  Preequillibrium config. Default set to -1
   * @param initialNeutronHoleNumber Preequillibrium config. Default set to -1
   * @param initialProtonNumber Preequillibrium config. Default set to -1
   * @param initialProtonHoleNumber Preequillibrium config. Default set to -1
   */
  DecayController(int Z, int A, double jInitial, int piInitial, double energy,
		 int initialNeutronNumber = -1, int initialNeutronHoleNumber = -1,
		 int initialProtonNumber = -1, int initialProtonHoleNumber = -1) :
    Z_(Z), A_(A), piInitial_(piInitial), initialNeutronNumber_(initialNeutronNumber), initialNeutronHoleNumber_(initialNeutronHoleNumber),
    initialProtonNumber_(initialProtonNumber), initialProtonHoleNumber_(initialProtonHoleNumber),
    jInitial_(jInitial), energy_(energy) {};

  bool Decay(double&,double&,double&,double&,double&,double&,double&,double&);

  std::vector<DecayProduct> DecayProducts() const { return decayProducts_; };

  double Energy(){return energy_;}  /**< Getter for energy_*/
 
  /**
   * @brief Print decay informations to std::cout
   */
  void PrintDecays();

 private:
 
  /**
   * @brief 
   * @param daughterZ Charge of the daughter
   * @param daughterA Mass of the daughter
   * @param parentZ Charge of the parent
   * @param parentA Mass of the parent
   * @param daughterEnergy Energy of the daughter
   * @param daughterJ Spin of the daughter
   * @param daughterPi Parity of the daughter
   * @param decayEnergy Decay energy \f$E_d\f$
   * @param parentBeta Velocity of the parent as a TVector3
   */
  void CalcKinematics(int daughterZ, int daughterA, int parentZ, int parentA, double daughterEnergy, double daughterJ, int daughterPi, double decayEnergy, TVector3 &parentBeta);
 
 private:
  
  static constexpr double thresholdEnergy_=0.05;
  int Z_;
  int A_;
  int piInitial_;
  int initialNeutronNumber_;
  int initialNeutronHoleNumber_;
  int initialProtonNumber_;
  int initialProtonHoleNumber_;
  double jInitial_;
  double energy_;
  std::vector<DecayProduct> decayProducts_;
};
