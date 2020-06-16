#pragma once

/**
 * @brief Basic class for transmission functions. Parent of ParticleTransmissionFunc and GammaTransmissionFunc.
 */
class TransmissionFunc {
 public:
 /**
  * @brief Constructor of a TransmissionFunc.
  * @param z2
  * @param m2
  * @param jInitial
  * @param piInitial
  * @param jFinal
  * @param piFinal
  * @param maxL
  * @param TWFC
  * @param uTWFC
  * @param uTWSFC
  * @param previous
  */
  TransmissionFunc (int z2, int m2, double jInitial, int piInitial,
		    double jFinal, int piFinal, double maxL,
		    double TWFC,
		    double uTWFC,
		    double uTWSFC,
		    TransmissionFunc* previous) :
    z2_(z2),m2_(m2),piInitial_(piInitial),piFinal_(piFinal),
    jInitial_(jInitial),jFinal_(jFinal), maxL_(maxL),
    TWFC_(TWFC),uTWFC_(uTWFC), uTWSFC_(uTWSFC),
    previous_(previous) {};  

  /**
   * @brief Destructor.
   */
  virtual ~TransmissionFunc() {};

  /**
   * @brief operator()
   */
  virtual double operator()(double) = 0;

  virtual bool IsValid() = 0;

 protected:
  int z2_;  /**< charge numbere*/
  int m2_; /**< mass number */
  int piInitial_;
  int piFinal_;
  double jInitial_;
  double jFinal_;
  double maxL_;   /**< maximum angular momentum*/
  double TWFC_; /**< Total width for correction*/
  double uTWFC_; /**< uncorrected total width for correction*/
  double uTWSFC_; /**< uncorrected total width squared for correction*/
  TransmissionFunc* previous_;  /**< previous transmission func; needed for width fluctuation correction*/
};


