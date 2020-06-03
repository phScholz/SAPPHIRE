/**
 * @file DecayProduct.h
 * @brief Declaration of DecayData and DecayProduct
 */
#pragma once

/**
 * @brief Class for the decay data
 */
class DecayData {
  double energy_;
  double spin_;
  int parity_;
  double neutronEntranceWidth_;
  double protonEntranceWidth_;
  double alphaEntranceWidth_;
  double gammaEntranceWidth_;
  double neutronTotalWidth_;
  double protonTotalWidth_;
  double alphaTotalWidth_;
  double gammaTotalWidth_;

 public:
  /**
   * @brief Simple DecayData class constructor. Initializes members with 0.
   */
  DecayData() : energy_(0.), spin_(0), parity_(0), neutronEntranceWidth_(0.), protonEntranceWidth_(0.),alphaEntranceWidth_(0.), gammaEntranceWidth_(0.), neutronTotalWidth_(0.), protonTotalWidth_(0.), alphaTotalWidth_(0.),gammaTotalWidth_(0.) {};
  
  /**
   * @brief DecayData class constructor. Initializes members with function parameters.
   * @param energy
   * @param spin
   * @param parity
   * @param neutronEntranceWidth
   * @param protonEntranceWidth
   * @param alphaEntranceWidth
   * @param gammaEntranceWidth
   * @param neutronTotalWidth
   * @param protonTotalWidth
   * @param alphaTotalWidth
   * @param gammaTotalWidth
   */
  DecayData(double energy, double spin, int parity, double neutronEntranceWidth,double protonEntranceWidth,
	   double alphaEntranceWidth, double gammaEntranceWidth,double neutronTotalWidth,
	   double protonTotalWidth, double alphaTotalWidth, double gammaTotalWidth) : energy_(energy), spin_(spin), parity_(parity), neutronEntranceWidth_(neutronEntranceWidth),
    protonEntranceWidth_(protonEntranceWidth), alphaEntranceWidth_(alphaEntranceWidth),
    gammaEntranceWidth_(gammaEntranceWidth), neutronTotalWidth_(neutronTotalWidth), protonTotalWidth_(protonTotalWidth),
    alphaTotalWidth_(alphaTotalWidth),gammaTotalWidth_(gammaTotalWidth) {};  

  double energy() const {return energy_;}; /**< Getter energy_*/
  double spin() const {return spin_;}; /**< Getter spin_*/
  int parity() const {return parity_;}; /**< Getter parity_*/
  
  double neutronEntranceWidth() const {return neutronEntranceWidth_;};/**< Getter neutronEntranceWidth_ */
  double protonEntranceWidth() const {return protonEntranceWidth_;};/**< Getter protonEntranceWidth_ */
  double alphaEntranceWidth() const {return alphaEntranceWidth_;};/**< Getter  alphaEntranceWidth_*/
  double gammaEntranceWidth() const {return gammaEntranceWidth_;};/**< Getter gammaEntranceWidth_*/
  double neutronTotalWidth() const {return neutronTotalWidth_;};/**< Getter neutronTotalWidth_ */
  double protonTotalWidth() const {return protonTotalWidth_;};/**< Getter protonTotalWidth_*/
  double alphaTotalWidth() const {return alphaTotalWidth_;};/**< Getter alphaTotalWidth_*/
  double gammaTotalWidth() const {return gammaTotalWidth_;};/**< Getter gammaTotalWidth_*/
};

/**
 * @brief Class for Decay Products
 */
class DecayProduct {
 
 public:
  /**
   * @brief Simple Constructor for DecayProduct class object.
   */
  DecayProduct() :
   Z_(0), A_(0), J_(0.), Pi_(0), excitationEnergy_(0.),
   fragmentEnergyCM_(0.), fragmentEnergy_(0.),
   fragmentMomentumX_(0.),fragmentMomentumY_(0.),
   fragmentMomentumZ_(0.), particleType_(0), 
   particleThetaCM_(0.), particlePhiCM_(0.),
   particleEnergyCM_(0.), particleEnergy_(0.), 
   particleMomentumX_(0.), particleMomentumY_(0.),
   particleMomentumZ_(0.) {};
  
  /**
   * @brief DecayProduct constructor 
   */
  DecayProduct(int Z, int A, double J, int Pi, 
	       double excitationEnergy, double fragmentEnergyCM,
	       double fragmentEnergy, double fragmentMomentumX,
	       double fragmentMomentumY, double fragmentMomentumZ,
	       int decayType, double particleThetaCM , double particlePhiCM, 
	       double particleEnergyCM, double particleEnergy, 
	       double particleMomentumX, double particleMomentumY,  
	       double particleMomentumZ) :
    Z_(Z), A_(A), J_(J), Pi_(Pi),excitationEnergy_(excitationEnergy),
    fragmentEnergyCM_(fragmentEnergyCM), fragmentEnergy_(fragmentEnergy),
    fragmentMomentumX_(fragmentMomentumX), fragmentMomentumY_(fragmentMomentumY),
    fragmentMomentumZ_(fragmentMomentumZ), particleType_(decayType), 
    particleThetaCM_(particleThetaCM), particlePhiCM_(particlePhiCM),
    particleEnergyCM_(particleEnergyCM), particleEnergy_(particleEnergy),
    particleMomentumX_(particleMomentumX),particleMomentumY_(particleMomentumY),
    particleMomentumZ_(particleMomentumZ) {};
    
  int Z_;
  int A_;
  int Pi_;
  int particleType_;
  double J_;
  double excitationEnergy_;
  double fragmentEnergyCM_;
  double fragmentEnergy_;
  double fragmentMomentumX_;
  double fragmentMomentumY_;
  double fragmentMomentumZ_;
  double particleThetaCM_;
  double particlePhiCM_;
  double particleEnergyCM_;
  double particleEnergy_;
  double particleMomentumX_;
  double particleMomentumY_;
  double particleMomentumZ_;
};

class DecayEvent{
  public:
    DecayData  Data;
    std::vector<DecayProduct> Products;
};