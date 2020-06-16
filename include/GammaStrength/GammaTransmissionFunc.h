/**
 * @file GammaTransmissionFunc.h
 * @date 2020-04-27
 * @brief Declaration of the base class for Gamma Strength Functions
 */
#pragma once
#include "../TransmissionFunc.h"
#include "Constants.h"
#include "Databases/NuclearMass.h"
#include "LevelDensity/LevelDensityFormula.h"
#include <math.h>

/**
 * @brief Class to store GDR parameters for an isotope
 * @details 
 * The values and errors of giant dipole resonance (GDR) parameters are presented which were obtained by a fit of the theoretical photoabsorption cross sections to the experimental data for 121 nuclides from 12-C through 239-Pu. The values and errors of the shape parameters of the Lorentzian-like curves corresponding to the giant dipole resonance excitation are presented.[1-8]
 *
 *  References:
 * 
 * - [1] J. Kopecky, in Handbook for calculations of nuclear reaction data. Reference Input Parameter Library (RIPL), IAEA-TEDOC-1034, 1998, Ch.6
 * - [2] T. Belgya, O. Bersillon, R. Capote, T. Fukahori, G. Zhigang, S. Goriely, M. Herman, A.V. Ignatyuk, S. Kailas. A. Koning, P. Oblozinsky, V. Plujko, P. Young. Handbook for calculations of nuclear reaction data: Reference Input Parameter Library-2, IAEA-TECDOC-1506, Vienna, 2006, Ch.7.
 * - [3] V.A. Plujko, I.M. Kadenko, E.V. Kulich, S. Goriely, O.I. Davidovskaya, O.M. Gorbachenko, in Proceeding of Workshop on Photon Strength Functions and Related Topics, Prague, Czech Republic, June 17-20, 2007, Proceedings of Science, PSF07, 2008
 * - [4] S.S.Dietrich, B.L.Berman; At. Data Nucl. Data Tables., 199, 38(1988).
 * - [5] M.B. Chadwick, P. Oblozinsky, P.E. Hodgson, G. Reffo. Phys.Rev. C44(1991)814.
 * - [6] M.B.Chadwick, P.Oblozinsky, A.I.Blokhin, T.Fukahori, Y.Han, Y. O.Lee, M.N.Martins, S.F.Mughabghab, V.V.Varlamov, B.Yu, J.Zhang. Handbook on photonuclear data for application. Cross sections and spectra. IAEA TECDOC 1178, Vienna, 2000
 * - [7] Experimental Nuclear Reaction Data Library EXFOR
 * - [8] CERN Program Library, MINUIT (D506), Function Minimization and Error Analysis
 * 
 * See the [RIPL-3 README](https://www-nds.iaea.org/RIPL-3/gamma/gdr-parameters&errors-exp.readme)
 * for details.
 */
class GDRParameters {
 public:
 /**
  * @brief constructor; initializing members to zero.
  */
  GDRParameters() : eta_(0.) {
    E_[0]=0.;
    E_[1]=0.;
    W_[0]=0.;
    W_[1]=0.;
  };

  /**
   * @brief Constructor; Initializing members from parameters.
   */
 GDRParameters(double eta, double E1, double W1, double kSigmaGamma1,
	       double E2, double W2, double kSigmaGamma2) : 
   eta_(eta) {
    E_[0]=E1;
    E_[1]=E2;
    W_[0]=W1;
    W_[1]=W2;
    kSigmaGamma_[0] = kSigmaGamma1;
    kSigmaGamma_[1] = kSigmaGamma2;
  };
  double eta_; 
  double E_[2]; /**< energy positions of the peaks*/
  double W_[2]; /**< widths of the peaks*/
  double kSigmaGamma_[2]; /**< strength of the peaks*/
};


typedef std::tr1::unordered_map<MassKey, GDRParameters > GDRTable; /**< Table to store GDR parameters*/

/**
 * @brief Base class for gamma-ray strength functions. Child of TransmissionFunc
 */
class GammaTransmissionFunc : public TransmissionFunc {
 public:
    /**
     * @brief constructor of GammaTransmissionFunc
     * @param z2 Charge number of compound nucleus
     * @param m2 Mass number of compound nucleus
     * @param jInitial Initial spin \f$ J_i \f$
     * @param piInitial Initial parity \f$ \pi_i \f$
     * @param jFinal Final spin \f$ J_f \f$
     * @param piFInal Final parity \f$ \pi_f \f$
     * @param maxL Maximum difference in angular momentum \f$ dL_{max} \$f
     * @param TWFC 
     * @param uTWFC
     * @param uTWSFC
     * @param previous
     */
    GammaTransmissionFunc(int z2,int m2, double jInitial, int piInitial, double jFinal,
    int piFinal, double maxL, double TWFC, double uTWFC, double uTWSFC, TransmissionFunc* previous); 
    
    virtual ~GammaTransmissionFunc(){};
    
    /**
     * @brief Chooses the Gammastrength model depending on input parameters and maximum angular momentum difference.
     */
    static GammaTransmissionFunc* CreateGammaTransmissionFunc(int z2,int m2, double jInitial, int piInitial, double jFinal,
    int piFinal, double maxL, double TWFC, double uTWFC, double uTWSFC, TransmissionFunc* previous, double compoundE);
    
    bool IsValid() {
      return true;
    };

  double operator()(double);

  /**
   * @brief Virtual function to calculate strength function values.
   * @param energy Energy of the transitions.
   * @details
   * This function need to be defined for any strength function model.
   */
  virtual double CalcStrengthFunction(double energy) = 0;

  /**
   * @brief Function to read the GDR parameters from file. 
   * @param file Path to the file containing the GDR parameters.
   * @details
   * Will be called in the Initialize() function from Setup.cpp 
   * and will fill the GDRTable gdrParameters_.
   */
  static void InitializeGDRParameters(std::string file);

  /**
   * @brief Static function to set the EGDRType (E1)
   * @param type Model for the strength function
   */
  static void SetEGDRType(int type){egdrType_=type;};

  /**
   * @brief Static function to set the MGDRType (M1)
   * @param type Model for the strength function
   */
  static void SetMGDRType(int type){mgdrType_ = type;};

  /**
   * @brief Static function to set the EGQRType (E2)
   * @param type Model for the strength function
   */
  static void SetEGQRType(int type){egqrType_ = type;};

  /**
   * @brief Static function to toggle Porter-Thomas fluctuations.
   * @param toggle True or false.
   */
  static void SetPorterThomas(bool toggle){porterThomas_=toggle;};

  /**
   * @brief Set the gnorm parameter
   * @param x The value the transmission coefficient for gammas should be multiplied with*/
  static void SetGnorm(double x){gnorm_=x;}

 protected:
  static GDRTable gdrTable_; /**< Map to store GDR parameters*/
  static int egdrType_; /**< Type of E1 strength*/
  static int mgdrType_; /**< Type of M1 strength*/
  static int egqrType_; /**< Type of E2 strength*/
  static bool porterThomas_; /**< Toggle Porter Thomas fluctuations.*/
  GDRParameters gdrParameters_; /**< Current GDR parameters for the nucleus*/
  static double gnorm_; /**< total normalization of transmission coefficient for gammas*/
};

