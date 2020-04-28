/**
 * @file ParticleTransmissionFunc.h
 * @brief Declaration of ParticleTransmissionFuncs.
 */
#pragma once

#include "TransmissionFunc.h"
#include "NuclearMass.h"
#include "Constants.h"
#include "SLPair.h"
#include <map>


/**
 * @brief Class for the particle transmission function. Child of TransmissionFunc.
 */
class ParticleTransmissionFunc : public TransmissionFunc {
 public:

 /**
  * @brief Constructor of ParticleTransmissionFunc
  * @param z1 Charge of the projectile
  * @param m1 Mass of the projectile
  * @param z2 Charge of the target nucleus
  * @param m2 Mass of the target nucleus
  * @param jInitial Initial spin of the target nucleus
  * @param piInitial Initial parity of the target nucleus
  * @param jFinal Final spin of the nucleus
  * @param piFinal Final parity of the nucleus
  * @param maxL Maximum of considered angular momentums
  * @param totalWidthForCorrection
  * @param uncorrTotalWidthForCorrection
  * @param uncorrTotalWidthSqrdForCorrection
  * @param previous
  * @details
  * 1. Initializes member variables 
  * 2. Calculates reduced mass, but only if the masses can be found via NuclearMass:FindMass()
  * 3. Selects type of projectile from z1 and m1
  * 4. Return transmissionFunc depending the chosen input for formalism
  */
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
    
    if(!NuclearMass::FindMass(z1_,m1_,mass1) || !NuclearMass::FindMass(z2_,m2_,mass2)) redmass_=0.;
    else redmass_=mass1*mass2/(mass1+mass2)/uconv;
    
    if(z1_==0&&m1_==1) pType_=0;        //neutron
    else if(z1_==1&&m1_==1) pType_=1;   //proton
    else if(z1_==2&&m1_==4) pType_=2;   //alpha
    else pType_=-1;
  };

  /**
   * @brief Destructor of ParticleTransmissionFunc
   */
  virtual ~ParticleTransmissionFunc() {};

  /**
   * @brief Function to check whether the input is valid or not.
   * @details
   * Returns True only if the reduced mass is non-zero and the projectile a neutron, proton or alpha.
   */
  bool IsValid() {
    if(redmass_==0.||pType_==-1) return false;
    else return true;
  };

  /**
   * @brief Operator(double)
   * @param energy 
   * @returns 
   */
  double operator()(double energy);

  /**
   * @brief Operator(double,int)
   * @param energy
   * @param which
   */
  double operator()(double energy,int which);
  void CalcSLDependentFunctions(double,std::map<SLPair,double>&);

  /**
   * @brief Create a ParticleTransmissionFunc. Will be called from TransitionRateFunc().
   * @param z1 Charge number of isotope 1
   * @param m1 Mass number of isotope 1
   * @param z2 Charge number of isotope 2
   * @param m2 Mass number of isotope 2
   * @param jInitial Initial spin
   * @param piInitial Initial parity
   * @param jFinal Final spin
   * @param piFinal Final parity
   * @param spin,Spin of ejectile
   * @param parity Parity of ejectile
   * @param maxL maximum angular momentum
   * @param totalWidthForCorrection
   * @param uncorrTotalWidthForCorrection
   * @param uncorrTotalWidthSqrdForCorrection
   * @param previous
   * @return ParticleTransmissionFunc object.
   * @details
   * 1. Check the projectile according to z1 and m1
   * 2. Creates a new TransmissionFunc object based on the input-parameter of the formalism
   * 3. Returns the created TransmissionFunc
   */
  static ParticleTransmissionFunc* CreateParticleTransmissionFunc(int z1,
  int m1, int z2, int m2, double jInitial, int piInitial, double jFinal, 
  int piFinal, double spin, int parity, double maxL, double totalWidthForCorrection,
  double uncorrTotalWidthForCorrection, double uncorrTotalWidthSqrdForCorrection,
  TransmissionFunc* previous);
  
  static void SetAlphaFormalism(int formalism) {alphaFormalism_ = formalism;};      /**< Setter for AlphaFormalism*/
  static void SetNeutronFormalism(int formalism) {neutronFormalism_ = formalism;};  /**< Setter for NeutronFormalism*/
  static void SetProtonFormalism(int formalism) {protonFormalism_ = formalism;};    /**< Setter for ProtonFormalism*/
  static void SetPorterThomas(bool yn){ porterThomas_=yn;} /**< Setter for porterThomas_*/
 
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


