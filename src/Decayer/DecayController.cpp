/**
 * @file DecayController.cpp
 * @brief Logic of the DecayController class
 */

#include "Decayer/DecayController.h"
#include "Constants.h"
#include "Databases/NuclearMass.h"
#include "Decayer/PreEqDecayer.h"
#include <iostream>
#include <iomanip>
#include <TVector3.h>
#include <stdlib.h>
#include <omp.h>

extern unsigned int randomSeed[12];


bool DecayController::Decay(double& neutronEntrance,
			                      double& protonEntrance,
			                      double& alphaEntrance,
			                      double& gammaEntrance,
			                      double& neutronTotalWidth,
			                      double& protonTotalWidth,
			                      double& alphaTotalWidth,
			                      double& gammaTotalWidth) 
  {
    bool returnValue = true;
    int parentA = A_;
    int parentZ = Z_;
    int parentPi = piInitial_;
    double parentJ = jInitial_;
    double parentEnergy = energy_;
    TVector3 parentBeta = TVector3(0.,0.,0.);  
    int parentNeutronNumber = initialNeutronNumber_;
    int parentNeutronHoleNumber = initialNeutronHoleNumber_;
    int parentProtonNumber = initialProtonNumber_;
    int parentProtonHoleNumber = initialProtonHoleNumber_;
    
    Decayer *initialEqDecayer = new Decayer(parentZ,parentA,parentJ,parentPi,parentEnergy);
    //initialEqDecayer -> PrintCDF();
    neutronEntrance=initialEqDecayer->NeutronEntranceWidth();
    protonEntrance=initialEqDecayer->ProtonEntranceWidth();
    alphaEntrance=initialEqDecayer->AlphaEntranceWidth();
    gammaEntrance=initialEqDecayer->GammaEntranceWidth();
    gammaTotalWidth=initialEqDecayer->GammaTotalWidth();
    neutronTotalWidth=initialEqDecayer->NeutronTotalWidth();
    alphaTotalWidth=initialEqDecayer->AlphaTotalWidth();
    protonTotalWidth=initialEqDecayer->ProtonTotalWidth();

    delete initialEqDecayer;

    bool decay = true;
    if(initialNeutronNumber_>0||initialNeutronHoleNumber_>0||initialProtonNumber_>0||initialProtonHoleNumber_>0)
    {
      bool preEq = true;
      int numSameExciton = 0;
      
      while(preEq) {
        PreEqDecayer* decayer = new PreEqDecayer(parentNeutronNumber,parentNeutronHoleNumber,
					                                       parentProtonNumber,parentProtonHoleNumber,
					                                      parentZ,parentA,parentJ,parentPi,parentEnergy);
        int daughterZ,daughterA,daughterPi;
        int daughterNeutronNumber;
        int daughterNeutronHoleNumber;
        int daughterProtonNumber;
        int daughterProtonHoleNumber;
        double daughterJ, daughterEnergy, decayEnergy;
        
        if(decayer->Decay(daughterZ,daughterA,daughterJ,daughterPi,daughterNeutronNumber,daughterNeutronHoleNumber,daughterProtonNumber,daughterProtonHoleNumber,daughterEnergy,decayEnergy))
        {
	        int daughterTotalExcitonNumber = daughterNeutronNumber+daughterNeutronHoleNumber+daughterProtonNumber+daughterProtonHoleNumber;
          int parentTotalExcitonNumber = parentNeutronNumber+parentNeutronHoleNumber+parentProtonNumber+parentProtonHoleNumber;

          //decayer->PrintCDF();
	        /*preEq = false;
	        decay = false;
	        continue;*/

	        if(daughterTotalExcitonNumber==parentTotalExcitonNumber) numSameExciton++;
	        else numSameExciton=0;
	
          if(numSameExciton>2||(daughterNeutronNumber>5&&daughterProtonNumber>5)|| daughterTotalExcitonNumber==1) 
	          preEq=false;

	        if(daughterEnergy<=thresholdEnergy_) {
	          decayEnergy+=daughterEnergy;
	          daughterEnergy = 0.;
	          preEq = false;
	          decay = false;
	        }

	        CalcKinematics(daughterZ,daughterA,parentZ,parentA,daughterEnergy,daughterJ,daughterPi,decayEnergy,parentBeta);

	        if(preEq) {
	          parentNeutronNumber=daughterNeutronNumber;
	          parentNeutronHoleNumber=daughterNeutronHoleNumber;
	          parentProtonNumber=daughterProtonNumber;
	          parentProtonHoleNumber=daughterProtonHoleNumber;
	        }

	        if(decay) {
	          parentZ = daughterZ;
	          parentA = daughterA;
	          parentJ = daughterJ;
	          parentPi = daughterPi;
	          parentEnergy = daughterEnergy;
	        }
        } else {    
	        returnValue = false;
	        decay = false;
	        preEq=false;
        }
        delete decayer;
      } 
    }

    while (decay) {
      Decayer *decayer = new Decayer(parentZ,parentA,parentJ,parentPi,parentEnergy);
      int daughterZ,daughterA,daughterPi;
      double daughterJ, daughterEnergy, decayEnergy;
    
      if(decayer->Decay(daughterZ,daughterA,daughterJ,daughterPi,daughterEnergy,decayEnergy))
      {
        
/*      decayer->PrintFunctions();*/
        //decayer->PrintCDF();
        /*decay = false;
        continue; */
                

        if(daughterEnergy<=thresholdEnergy_){
	        daughterEnergy = 0.;
	        decayEnergy+=daughterEnergy;
	        decay = false;
        }

        CalcKinematics(daughterZ,daughterA,parentZ,parentA,daughterEnergy,daughterJ,daughterPi,decayEnergy,parentBeta);

        if(decay){
	        parentZ = daughterZ;
	        parentA = daughterA;
	        parentJ = daughterJ;
	        parentPi = daughterPi;
	        parentEnergy = daughterEnergy;
        }
      } else {
        returnValue = false;
        decay = false;
      }
      delete decayer;
    }
    
    return returnValue;
  }

void DecayController::PrintDecays() {
  /**
   * 1. Print the header.
   */
  std::cout << std::endl
	    << std::setw(5)  << "Z"  << ' '
	    << std::setw(5)  << "A"  << ' '
	    << std::setw(5)  << "J"  << ' '
	    << std::setw(5)  << "Pi" << ' '
	    << std::setw(15) << "Heavy Exc. En." << ' '
	    << std::setw(5)  << "Type" << ' '
	    << std::setw(14) << "Light En. CM" << ' '
	    << std::setw(14) << "Heavy En. CM" << ' '
	    << std::setw(14) << "Light En. Lab" << ' '
	    << std::setw(14) << "Heavy En. Lab" << ' '
	    << std::setw(0)  << std::endl;

  /**
   * 2. For each element in decayProducts, print Z, A, J, Pi, excitation energy, particle Type
   * and the center-of-mass and laboratory energies of the particle and the fragment.
   */
  for(int i = 0; i<decayProducts_.size();i++)
  {
    std::cout << std::setw(5)  << decayProducts_[i].Z_ << ' '
	            << std::setw(5)  << decayProducts_[i].A_ << ' '
	            << std::setw(5)  << decayProducts_[i].J_ << ' '
	            << std::setw(5)  << decayProducts_[i].Pi_ << ' '
	            << std::setw(15) << decayProducts_[i].excitationEnergy_ << ' ';
   
    if(decayProducts_[i].particleType_==0)
    { 
      std::cout << std::setw(5)  << 'g' << ' ';
    }
    else if(decayProducts_[i].particleType_==1)
    {
	    std::cout << std::setw(5)  << 'n' << ' ';
    }
    else if(decayProducts_[i].particleType_==2)
    {
	    std::cout << std::setw(5)  << 'p' << ' ';
    }
    else if(decayProducts_[i].particleType_==3)
    {
	    std::cout << std::setw(5)  << 'a' << ' ';      
    }
    else std::cout << std::setw(5)  << '?' << ' ';

    std::cout << std::setw(14) << decayProducts_[i].particleEnergyCM_ << ' '
	            << std::setw(14) << decayProducts_[i].fragmentEnergyCM_ << ' '
	            << std::setw(14) << decayProducts_[i].particleEnergy_ << ' '
	            << std::setw(14) << decayProducts_[i].fragmentEnergy_ << ' '
              << std::setw(0)  << std::endl;
  }
}

void DecayController::CalcKinematics(int daughterZ, int daughterA, int parentZ, int parentA, double daughterEnergy, double daughterJ, int daughterPi, double decayEnergy, TVector3 &parentBeta) {
  /**
   * 1. Define the decay type, from parent and daughter mass and charge.
   *  - Type 0 : Gamma-decay
   *  - Type 1 : Neutron-decay
   *  - Type 2 : Proton-decay
   *  - Type 3 : Alpha-decay
   *  - Type 4 : Beta- decay (Placeholder for now)
   *  - Type 5 : Beta+ decay (Placeholder for now)
   */
  int decayType;
  if((parentA-daughterA)==0&&(parentZ-daughterZ)==0) decayType=0;
  else if((parentA-daughterA)==1&&(parentZ-daughterZ)==0) decayType=1;
  else if((parentA-daughterA)==1&&(parentZ-daughterZ)==1) decayType=2;
  else if((parentA-daughterA)==4&&(parentZ-daughterZ)==2) decayType=3;
  else if((parentA-daughterA)==0&&(parentZ-daughterZ)==1) decayType=4;
  else if((parentA-daughterA)==0&&(parentZ-daughterZ)==-1) decayType=5;

  /** 
   * 2. Looking for the masses of the daughter and the ejectile in the MassTable. 
   * Mass of daughter is increased by daughterEnergy.
   */
  double m1,m2;
  if(!NuclearMass::FindMass(daughterZ,daughterA,m1)) {
    std::cout << "Unknown masses requested.  Aborting." << std::endl;
    exit(1);
  }
  if(decayType==0) m2=0.;
  else if(!NuclearMass::FindMass(parentZ-daughterZ,parentA-daughterA,m2)) {
    std::cout << "Unknown masses requested.  Aborting." << std::endl;
    exit(1);
  }
  m1+=daughterEnergy;

  /**
   * 3. Calculate total mass \f[ m_{Total}=m_1+m_2 \f] and total center of mass momentum 
   * \f[ p_{cm} = \frac{1}{1+E_d/m_{Total}} \cdot \sqrt{ E_d^4 / (4 m_{Total}^2) + E_d^3/m_{Total} + E_d^2 (1+\frac{m_1 m_2}{m_{Total}^2}) + 2 E_d (\frac{m_1 m_2}{m_{Total}})} \f]
   */
  double mTotal = m1+m2;
  double momentumCM = 1./(1.+decayEnergy/mTotal)* sqrt(pow(decayEnergy,4.)/4./pow(mTotal,2.)+ pow(decayEnergy,3.)/mTotal+ pow(decayEnergy,2.)*(1.+m1*m2/pow(mTotal,2.))+ 2.*decayEnergy*m1*m2/mTotal);

  /**
   * 4. \f$\phi_{cm} \f$ and \f$ \theta_{cm} \f$ are drawn randomly and center-of-mass vectors are created for the 
   * different momentums:
   *    - decay direction \f[ \vec{e}_d = \begin{pmatrix} \sin{\theta_{cm}}\cos{\pi_{cm}} \\\ \sin{\theta_{cm}}\sin{\pi_{cm}} \\\ \cos{\theta_{cm}} \end{pmatrix}  \f]
   *    - total center of mass decay momentum: \f[ p_{cm}^{decay} =  p_{cm} \cdot \vec{e}_d  \f]
   *    - total center of mass daughter momentum: \f[ p_{cm}^{daughter} = - p_{cm} \cdot \vec{e}_d  \f]
   *    - total center of mass decay energy: \f[ E_{cm}^{decay} = \sqrt{ p_{cm}^2 + m_2^2}  \f]
   *    - total center of mass daughter energy: \f[ E_{cm}^{daughter} = \sqrt{ p_{cm}^2 + m_1^2}  \f]
   *    - The decay laboratory momentum is calculated from the parentBeta and the total center-of-mass decay energy
   *    - The daughter laboratory momenutm is calculated from the parentBeta and the daughter center-of-mass energy and momentum
   *    - The decay and daughter energies are transformed to the laboratory frame.
   */
  double phiCM = 2.*pi*double(rand_r(&randomSeed[omp_get_thread_num()]))/double(RAND_MAX);
  double thetaCM = acos(2.*double(rand_r(&randomSeed[omp_get_thread_num()]))/double(RAND_MAX)-1.);

  TVector3 decayDirectionCM = TVector3(sin(thetaCM)*cos(phiCM),sin(thetaCM)*sin(phiCM),cos(thetaCM));
  TVector3 decayMomentumCM = momentumCM*decayDirectionCM;
  TVector3 daughterMomentumCM = -momentumCM*decayDirectionCM;

  double decayEnergyTotalCM = sqrt(momentumCM*momentumCM+m2*m2);
  double daughterEnergyTotalCM = sqrt(momentumCM*momentumCM+m1*m1);

  double parentGamma = 1./sqrt(1.-parentBeta.Mag2());

  TVector3 decayMomentum = (parentBeta.Mag2()==0.) ? decayMomentumCM : 
    decayMomentumCM + parentBeta*(parentGamma*decayEnergyTotalCM+
				  (parentGamma-1.)*parentBeta.Dot(decayMomentumCM)/
				  parentBeta.Mag2());

  double transformDecayKineticEnergy = parentGamma*(decayEnergyTotalCM+
							parentBeta.Dot(decayMomentumCM))-m2;

  TVector3 daughterMomentum = (parentBeta.Mag2()==0.) ? daughterMomentumCM : 
    daughterMomentumCM + parentBeta*(parentGamma*daughterEnergyTotalCM+
					 (parentGamma-1.)*parentBeta.Dot(daughterMomentumCM)/
					 parentBeta.Mag2());

  double transformDaughterKineticEnergy = parentGamma*(daughterEnergyTotalCM+
						       parentBeta.Dot(daughterMomentumCM))-m1;

  /**
   * 5. A DecayProduct is generated and pushed back to the member decayProducts_
   */
  decayProducts_.push_back(DecayProduct(daughterZ,daughterA,daughterJ,daughterPi,
					daughterEnergy,daughterEnergyTotalCM-m1,
					transformDaughterKineticEnergy,
					daughterMomentum.X(),daughterMomentum.Y(),daughterMomentum.Z(),
					decayType,thetaCM,phiCM,
					decayEnergyTotalCM-m2,transformDecayKineticEnergy,
					decayMomentum.X(),decayMomentum.Y(),decayMomentum.Z()));

  parentBeta=1./(transformDaughterKineticEnergy+m1)*daughterMomentum;
}
