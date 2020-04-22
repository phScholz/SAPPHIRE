#pragma once
#include <vector>
#include <map>
#include <string>
#include "SapphireInput.h"
#include "SpinRatePair.h"
#include "Decayer.h"

typedef std::vector<std::pair<Decayer*,std::vector<SpinRatePair*> > > DecayerVector;

typedef std::pair<int,double> int_double_pair;

struct int_double_pair_compare {
  bool operator()(const int_double_pair &lhs, int_double_pair const  &rhs) {
        return lhs.second < rhs.second;
  }
};


/**
 * @brief A class which represents an object of cross section data for all reactions on a nucleus at one energy.
 */
class CrossSectionValues {
 public:
 CrossSectionValues(double gamma, double neutron, double proton, double alpha,
		    double gammaStellar, double neutronStellar, double protonStellar, double alphaStellar) :
  gamma_(gamma),neutron_(neutron),proton_(proton),alpha_(alpha),
  gammaStellar_(gammaStellar),neutronStellar_(neutronStellar),protonStellar_(protonStellar),alphaStellar_(alphaStellar) {};
  double gamma_;
  double neutron_;
  double proton_;
  double alpha_;
  double gammaStellar_;
  double neutronStellar_;
  double protonStellar_;
  double alphaStellar_;
};

/**
 * @brief Class to perform the calculations for cross sections and reaction rates
 */
class CrossSection {
 public:
  /**
  * @brief The original constructor of the CrossSection class.
  * @param Z  Atomic number of the target nucleus.
  * @param A  Mass number of the target nucleus.
  * @param pType  Integer indicating the projectile type
  * @param energyFile The path to a file containing energies, or ""
  * @param forRates A boolean indicating whether rates should be calculated or not.
  * @param entranceState  The number of the initial state, e.g. 0 = groundstate
  * @param exitStates Std::vector<int> which contains the number of levels for which partial cross section should be calculated.
  * @details The constructor is responsible to setting up all internal variables needed for the calculation of cross sections.
  * The following steps will be taken:
  * 1. FindInitialState() will be called, which finds the entrance state in the list of known levels and stores the properties of 
  * spin and parity of this state into private member variables (CrossSection::groundstateJ_ and groundstatePi_).  
  * 2. PreSetCompound() is called to set the internal mass and atomic numbers for the compound nucleus and to calculate the prefactor of the cross section.
  * 3. Then the method InitializeSeperationEnergies() will be called which initializes the seperationEnergy_ as well as the 
  * properties of the exitStates.
  * 4. The validity of the exitStates will be checked via CheckChannels().
  * 5. If the input parameter was set to calculate reaction rates, then the normalized partition function \f$ G_0(T) = \frac{G(T)}{g_0} \f$ will be calculated via CalcPartitionFunc().
  * 6. CalcAllowedJPi() is called to fill the std::vector `allowedJPi_` with all allowed combinations of spin and parity for compoundstates accessible via the selection rules. Additionally, maps for the transmission coefficients are initialized.
  */
  CrossSection(int Z,int A,int pType,std::string energyFile,bool forRates,int entranceState = 0, std::vector<int> exitStates = std::vector<int>(4,-1));

  /**
  * @brief Constructor on the basis of a SapphireInput object.
  * @param input  SapphireInput object which needs to be initialized in Module_CrossSection
  * @details See CrossSection().
  */
  CrossSection(SapphireInput & input);

  bool IsValid() const {return isValid_;}                                       /**< Getter for isValid_*/

  void Calculate();
  void PrintCrossSections();
  void PrintTransmissionTerms();
  std::pair<double,double> CalcAverageSWaveResWidth();
  std::pair<double,double> CalcAveragePWaveResWidth();
  std::pair<double,double> CalcAverageDWaveResWidth();
  void CalculateReactionRates(bool);
  void PrintReactionRates(bool);
  static void SetResidualGamma(bool residual) {residualGamma_=residual;};         /**< Setter for residualGamma_*/
  static void SetResidualNeutron(bool residual) {residualNeutron_=residual;};     /**< Setter for residualNeutron_*/
  static void SetResidualProton(bool residual) {residualProton_=residual;};       /**< Setter for residualProton_*/
  static void SetResidualAlpha(bool residual) {residualAlpha_=residual;};         /**< Setter for residualAlpha_*/
  static void SetCalculateGammaCutoff(bool calc) {calculateGammaCutoff_=calc;};   /**< Setter for calculateGammaCutoff_*/
  static void CreateTempVector();
  static void CreateMACSEnergiesVector();

 private:
  /**
  * @brief Read in the energyFile and push back the energies to the CrossSection::crossSections_.
  * @param energyFile std::string containing the path to the energyFile
  * @return true or false depending, if the energies can be read or not
  */
  bool FillEnergies(std::string energyFile);

  /**
  * @brief Method to fill the `std::vector<std::pair<double,int> > allowedJPi_`.
  * @param calcRates True or false depending, if rates should be calculated or not.
  * @return True if the size of allowedJPi_ is larger than 0; False otherwise.
  * @details Depending on spin selection rules, the vector `allowedJP_` is filled under consideration of the target nucleus and the projectile. Additionally, the vectors for the maps: 
  * - `entranceTrans_`
  * - `gExitTrans_`
  * - `nExitTrans_`
  * - `pExitTrans_`
  * - `aExitTrans_`
  * are initialized using the size of `allowedJPi_`
  */
  bool CalcAllowedJPi(bool calcRates);
  bool CalcDecayerVector(double,DecayerVector&,bool forAverageWidth=false);
  
  /**
  * @brief Method to quickly check if the decay into the given exitChannels is energetically allowed.
  */
  void CheckChannels();
  
  /**
  * @brief Method to find the attributes of the entrance state and pass them to private member variables by reference.
  * @details The known levels of the target nucleus are obtained from NuclearLeves::FindLevels() using CrossSection::Z_ and CrossSection::A_ 
  * and stored into a std::vector<Level> knonLevels.
  * If the size of the known levels is larger than the number of CrossSection::entranceState_, then the 
  * properties (spin and parity) are stored in CrossSection::groundStateJ_ and CrossSection::groundStatePi_.
  * If the initial state cannot be found in the known levelscheme, the function returns false, otherwise true.
  * 
  * @return True or false depending on if the initial state is known or not.
  */
  bool FindInitialState();

  /**
  * @brief Method to find the compound mass, charge and preFactor and writes it into the member variables.
  * @details This method is responsible to setting the mass and charge number of the compund nucleus for different projectiles. 
  * Moreover, it calculates the preFactor_ in barn for the cross section calculation which is given by
  * \f[ c= \frac{(\hbar c)^2}{100} \cdot \frac{\pi}{2} \cdot \frac{A+1}{A \cdot u} \cdot \frac{1}{2\cdot(2J+1)} \f]
  * Here is \f$ A \f$ and \f$ J\f$ the mass number and spin of the target nucleus, respectively. The factor 100 is for the conversion of fm² to barn. 
  * @return true or false, depending if the given reaction qvalue is valid.
  */
  bool PreSetCompound();

  void InitializeSeperationEnergies();

  void CalculateEnergyGrid();

  /**
  * @brief Calculation of the normalized partition function of the nucleus as a function of temperatur.
  * The goal of this method is to calculate the partition function of the nucleus as a function of temperature.
  * The temperatures are obtained by the two std::vectors `rateTemps_` and `macsEnergies_`.
  * The normalized partition function is defined as 
  * \f[ G_0 = \frac{G(T)}{g_0} = \frac{1}{g_0} \left[ \sum_{\mu} g_{\mu} e^{-\frac{E^x_{\mu}}{kT}} + \int_{E^x_{\mu last}}^\infty
  * \sum_{J, \pi} g_J e^{-\frac{\epsilon}{kT} \rho(\epsilon,J,\pi)} d\epsilon \right]. \f] 
  * The first part is a simple sum over the boltzmann factors for degenerated states of the known levels of the nucleus.
  * The second part is the continuum extension using an integration over a level density. At the moment the RauscherLevelDensity is used for this.
  *
  * More information can be found in, e.g.,:
  * - T. Rauscher, *The path to improved reaction rates for astrophysics*, Internation Journal of Modern Physics E **20** (2011) 1071. DOI:[10.1142/S021830131101840X](https://doi.org/10.1142/S021830131101840X).
  * - C. Iliadis, *Nuclear Physics of Stars*, WILEY-VCG Verlag GmbH & Co. KGaA, 2007, ISBN: 978-3-527-40602-9. 
  *
  *For the numerical integration the GNU Scientific Library is used. 
  * Information about the integration method can be found here: [Numerical Integration with GSL](https://www.gnu.org/software/gsl/doc/html/integration.html)
  */
  void CalcPartitionFunc();

 private:
  static bool residualGamma_;                 /**< Bool if residualGamma should be toggled*/
  static bool residualNeutron_;               /**< Bool if residualNeutron should be toggled*/
  static bool residualProton_;                /**< Bool if residualProton should be toggled*/
  static bool residualAlpha_;                 /**< Bool if residualAlpha should be toggled*/
  static bool calculateGammaCutoff_;          
  constexpr static double minEnergy_ = 0.001;
  constexpr static double maxEnergy_ = 15.0;
  constexpr static double dE_ = .1;
  bool gammaCutoffSet_;
  bool isValid_;                              /**< Boolean to control, whether the input is correct or not*/
  bool energiesGiven_;
  bool calcRates_;
  bool verbose_= false;                       /**< Bool to increase the output for debugging purposes*/
  int Z_;                                     /**< The atomic number of the target nucleus*/
  int A_;                                     /**< The mass number of the target nucleus*/
  int pType_;                                 /**< Integer which defines the type of the projectile*/
  int compoundA_;                             /**< The mass number of the compound nucleus*/
  int compoundZ_;                             /**< The charge number of the compound nucleus*/
  int groundStatePi_;                         /**< The parit of the groundState*/
  int entranceState_;                         /**< The number of the entranceState in the level scheme*/
  double preFactor_;  	                      /**< The prefactor for the crosssection without energy. See PreSetCompound().*/
  double groundStateJ_;
  double seperationEnergy_;
  double skipEnergy_;
  double qValue_;
  std::vector<int> exitStates_;
  std::vector<double> specifiedExitSepE_;
  std::vector<double> specifiedExitJ_;
  std::vector<int> specifiedExitPi_;
  std::vector<bool> skipped_;
  std::vector<std::pair<double,int> > allowedJPi_;
  std::vector<std::pair<double,CrossSectionValues> > crossSections_;
  std::vector<std::pair<double,CrossSectionValues> > reactionRates_;
  std::map<int, std::vector<double> > entranceTrans_;
  std::map<int, std::vector<double> > gExitTrans_;
  std::map<int, std::vector<double> > nExitTrans_;
  std::map<int, std::vector<double> > pExitTrans_;
  std::map<int, std::vector<double> > aExitTrans_;
  std::vector<std::pair<double,double> > partFunc_;
  std::vector<std::pair<double,double> > partFuncMACS_;
  static std::vector<double> rateTemps_;
  static std::vector<double> macsEnergies_;
};

