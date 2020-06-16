/**
 * @file Decayer/Decayer.h
 * @brief Declaration of the Decayer class.
 */
#pragma once
#include <vector>
#include <cstdlib>
#include "SpinRatePair.h"
#include "TransitionRateFunc.h"
#include "SapphireInput.h"
#include "CDFEntry.h"
#include "Databases/NuclearLevels.h"

class Decayer{

  public:
  
    /**
     * @brief Constructor of the Decayer class
     * @param Z Charge/atomic number of the decaying nucleus
     * @param A Mass number of the decaying nucleus
     * @param jInitial The initial spin of the decaying nucleus
     * @param piInitial The initial parity of the decaying nucleus
     * @param energy The excitation energy of the decaying nucleus
     * @param totalWidthForCorrection By default 0.
     * @param uncorrTotalWidthForCorrection By default 0.
     * @param uncorrTotalWidthSqrdForCorrection By default 0.
     * @param widthCorrectedDecayer By default NULL.
     * @details 
     * 1. The following member variables will be initialized in the constructor:
     *  + Z_
     *  + A_
     *  + piInitial_
     *  + jInitial_
     *  + energy_
     *  + totalWidthForCorrection_
     *  + uncorrTotalWidthForCorrection_
     *  + uncorrTotalWidthSqrdForCorrection_
     *  + widthCorrectedDecayer_
     * 2. Calling of InitializeWidths() to initialize entrance and total widths for neutron, gamma, proton, and alpha are set to 0.
     * 3. Checking if the initial state is a bound state (BoundStateCheck()) of the nucleus and if yes -> build a cummulative probability distribution function (CDF) with BuildKnownCDF().
     * 4. If the state is not bound or no gamma-ray transitions are known from this state, proceed to use transitions from the continuum.
     * 5. The vec 
     */
    Decayer(int Z, int A, double jInitial, int piInitial, double energy, double totalWidthForCorrection = 0., double uncorrTotalWidthForCorrection = 0., double uncorrTotalWidthSqrdForCorrection = 0., Decayer* widthCorrectedDecayer=NULL);
    
    /**
     * @brief Constructor of Decayer on the basis of a SapphireInput onject
     * @param input A SapphireInput object with input-parameters
     * @details See Decayer()
     */
    Decayer(SapphireInput & input);

    /**
     * @brief Destructor of Decayer
     */
    ~Decayer();

    /**
     * @brief Drawing a randmomNumber and use that to determine the state to decay to via the CDF
     * @param Z Charge/atomic number of the daughter nucleus.
     * @param A Mass number of the daughter nucleus.
     * @param jFinal Spin of the daughter nucleus.
     * @param piFinal Parity of the daughter nucleus.
     * @param excitationEnergy Excitation energy of the daughter nucleus.
     * @param decayEnergy Decay energy.
     */
    bool Decay(int& Z,int& A,double& jFinal,int& piFinal,double& excitationEnergy,double& decayEnergy);
    

    /**
     * @brief Print levelDensity, transmissionFunc, transitionRateFunc, and totalLevelDensity to files.
     */
    void PrintFunctions();
    
    /**
     * @brief Writing the cummulative distribution function to cdf.out.
     * @details
     * The format of the file is
     * 
     *            Z A spin parity energy cdfvalue
     */
    void PrintCDF();

    void CorrectWidthFluctuations();

    /**
     * @brief Setter for the isCorssSection boolean
     */
    static void SetCrossSection(bool isCrossSection) {isCrossSection_=isCrossSection;};
    
    static void SetMaxL(double maxL) {maxL_=maxL;};       /**< Setter for maxL_*/
    static void SetAlphaDecay(bool x){alphaDecay_=x;}     /**<Setter for alphaDecay_*/
    static void SetProtonDecay(bool x){protonDecay_=x;}   /**<Setter for protonDecay_*/
    static void SetNeutronDecay(bool x){neutronDecay_=x;} /**<Setter for neutronDecay_*/
    static void SetGammaDecay(bool x){gammaDecay_=x;}     /**<Setter for gammaDecay_*/
    
    static double GetMaxL() {return maxL_;};            /**< Getter for maxL_*/
    
    friend class CrossSection;
    
    double NeutronEntranceWidth() const { return neutronEntrance_; } /**< Getter for neutronEntrance_*/
    double ProtonEntranceWidth() const { return protonEntrance_; }  /**< Getter for protonEntrance_*/
    double AlphaEntranceWidth() const { return alphaEntrance_; }    /**< Getter for alphaEntrance_*/
    double GammaEntranceWidth() const { return gammaEntrance_; }    /**< Getter for gammaEntrance_*/
    double GammaTotalWidth() const { return gammaTotalWidth_; }     /**< Getter for gammaTotalWidth_*/
    double NeutronTotalWidth() const { return neutronTotalWidth_; } /**< Getter for neutronTotalWidth_*/
    double AlphaTotalWidth() const { return alphaTotalWidth_; }     /**< Getter for alphaTotalWidth_*/
    double ProtonTotalWidth() const { return protonTotalWidth_; }   /**< Getter for protonTotalWidth_*/

  private:
    /**
     * @brief Initializes the member variables for partial and total widths.
     */
    inline void InitializeWidths();

    /**
     * @brief Looks up the QValues via NuclearMass::QValue() and writes it into qValueNeutron_, qValueProton_, and qValueAlpha_
     */
    inline void InitializeQValues();

    /**
     * 
     */
    void BuildCDF();
    
    /**
     * @brief Building the cummulative distribution function (CDF) for known levels of the decaying nucleus.
     * @param levelIndex The level number of a bound level until which the levelscheme should be built.
     * @param knownLevels A std::vector of Level entries of the known levels of the decaying nucleus
     * @return True or false depending on if gamma ray transitions can be found
     * @details
     * 1. For all available gamma-ray transitions of the known level, the probabilities will be summed up into totalIntegral_ as total decay width
     * 2. Values for the cummulative distribution function are calculated iteratively in a for loop over all gamma-ray transitions. 
     * The probability for the decay to the previous level is the offset of the probability for the decay to the next level.
     * That means, the probability to decay to the level with energy \f$ E_i \f$ is given by \f[ P(E_i) = CDF[E_i] - CDF[E_{i-1}] \f]
     * 3. For each gamma-ray a new CDFEntry is generated and pushed back to the cdf_ vector.
     */
    bool BuildKnownCDF(int levelIndex,std::vector<Level>& knownLevels);

    /**
     * @brief Check if the initial state is a bound state and if yes build the CDF via BuildKnownCDF().
     */
    inline bool BoundStateCheck();

    /**
     * @brief Adding SpinRatePairs entries for alpha decay.
     */
    inline void SetAlphaSpinRatePairs();

    /**
     * @brief Adding SpinRatePairs entries for proton or neutron decay.
     */
    inline void SetProtonNeutronSpinRatePairs();

    /**
     * @brief Adding SpinRatePairs entries for E1 or M1 transitions.
     */
    inline void SetE1M1SpinRatePairs();

    /**
     * @brief Adding SpinRatePairs entries for E2 transitions.
    */
    inline void SetE2SpinRatePairs();

   private:
    static bool isCrossSection_; /**< Boolean to indicate whether this Decayer is part of a CrossSection calculation or not.*/
    static double maxL_; /**< The maximum angular momentum*/
    static bool alphaDecay_; /**< Toggle alpha decay*/
    static bool protonDecay_; /**< Toggle proton decay*/
    static bool neutronDecay_; /**< Toggle neutron decay*/
    static bool gammaDecay_; /**< Toggle gamma decay*/
    int Z_; /**< The charge/atomic number of the decaying nucleus*/
    int A_; /**< The mass number of the decaying nucleus*/
    int piInitial_; /**< The parity of the intital state*/
    double jInitial_; /**< The spin of the initial state*/
    double energy_; /**< The excitation energy of the initial state*/
    double totalIntegral_; /**< The total decay width*/
    double totalIntegralSqrd_; /**< The square of the total decay widht*/
    double totalWidthForCorrection_;
    double uncorrTotalWidthForCorrection_;
    double uncorrTotalWidthSqrdForCorrection_;
    double neutronEntrance_; /**< neutron entrance width*/
    double protonEntrance_; /**< proton entrance width*/
    double gammaEntrance_;  /**< gamma entrance width*/
    double alphaEntrance_;  /**< alpha entrance width*/
    double neutronTotalWidth_;  /**< neutron total width*/
    double alphaTotalWidth_;  /**< alpha total width*/
    double protonTotalWidth_; /**< proton total width*/
    double gammaTotalWidth_;  /**< gamma total width*/
    double qValueProton_; /**< Qvalue for proton*/
    double qValueNeutron_;  /**< Qvalue for neutron*/
    double qValueAlpha_;  /**< Qvalue for alpha*/

    std::vector<SpinRatePair> spinRatePairs_; /**< Vector of Spin-Rate pairs*/
    std::vector<CDFEntry> cdf_; /**< Vector of CDFEntries from CDFEntry.h*/
    Decayer* widthCorrectedDecayer_; /**< Pointer to another Decayer for WFC*/
};
