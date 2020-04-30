/**
 * @file Setup.cpp
 * @brief Initialization of static members of the respective classes.
 * @details
 * In the original Sapphire code, input parameters were handled via static member
 * variables of the respective sub classes.
 * This is an easy way to set global parameters .. but maybe not the nicest.
 * However, it works like that. So, we continue to use it like that. Maybe once
 * in future, it would be better to template the respective classes for different models.
 */
#include "NuclearMass.h"
#include "GammaStrength/GammaTransmissionFunc.h"
#include "GammaStrength/D1MQRPA.h"
#include "NuclearLevels.h"
#include "Decayer.h"
#include "Sapphire_config.h"
#include "TransitionRateFunc.h"

#include "CrossSection.h"

#include "PreEqDecayer.h"
#include "ParticleTransmissionFunc.h"
#include "CoulFunc.h"
#include <iostream>
#include <gsl/gsl_errno.h>

#include "LevelDensity/LevelDensityTable.h"
#include "LevelDensity/LevelDensityHFB_BSk14.h"

//These are all static members of the respective classes

HFBTable LevelDensityHFB_BSk14::densityTable;
HFBCorrTable LevelDensityHFB_BSk14::corrTable;
double LevelDensityTable::c_;
double LevelDensityTable::d_;

QRPA_E1_Table D1MQRPA::e1Table;
QRPA_M1_Table D1MQRPA::m1Table;

bool CrossSection::residualGamma_;
bool CrossSection::residualNeutron_;
bool CrossSection::residualProton_;
bool CrossSection::residualAlpha_;
bool CrossSection::calculateGammaCutoff_;
std::vector<double> CrossSection::rateTemps_;
std::vector<double> CrossSection::macsEnergies_;

bool Decayer::alphaDecay_;
bool Decayer::protonDecay_;
bool Decayer::neutronDecay_;
bool Decayer::gammaDecay_;

bool Decayer::isCrossSection_;
bool PreEqDecayer::isCrossSection_;
double Decayer::maxL_;
double PreEqDecayer::maxL_;
double TransitionRateFunc::gammaCutoffEnergy_;
int TransitionRateFunc::nldmodel_;

ElementTable NuclearMass::elementTable_; 
MassTable NuclearMass::massTable_;
GDRTable GammaTransmissionFunc::gdrTable_;
LevelsTable NuclearLevels::levelsTable_;
int ParticleTransmissionFunc::alphaFormalism_;
int ParticleTransmissionFunc::neutronFormalism_;
int ParticleTransmissionFunc::protonFormalism_;
int GammaTransmissionFunc::egdrType_;
bool GammaTransmissionFunc::porterThomas_;
bool ParticleTransmissionFunc::porterThomas_;

/** 
 * @brief First Function which is called in the Sapphire main() in Sapphire.cpp. Setting default parameters
 * 
 */
void Initialize() {

  CrossSection::SetResidualGamma(true);
  CrossSection::SetResidualNeutron(false);
  CrossSection::SetResidualProton(false);
  CrossSection::SetResidualAlpha(false);
  CrossSection::SetCalculateGammaCutoff(true);
  CrossSection::CreateTempVector();
  CrossSection::CreateMACSEnergiesVector();
  
  //By default all decay channels are allowed if energetically possible
  Decayer::SetAlphaDecay(true);
  Decayer::SetProtonDecay(true);
  Decayer::SetNeutronDecay(true);
  Decayer::SetGammaDecay(true);

  Decayer::SetCrossSection(false);
  PreEqDecayer::SetCrossSection(false);
  
  Decayer::SetMaxL(8.);
  PreEqDecayer::SetMaxL(8.);
  TransitionRateFunc::SetGammaCutoffEnergy(10000.);
  TransitionRateFunc::NLDmodel(0);
  NuclearMass::InitializeElements();
  NuclearMass::InitializeMasses(sourceDirectory()+"/tables/masses/masses.dat");
  GammaTransmissionFunc::InitializeGDRParameters(sourceDirectory()+"/tables/gamma/ripl3_gdr_parameters.dat");
  NuclearLevels::InitializeLevels(sourceDirectory()+"/tables/levels/", sourceDirectory()+"/tables/spinod.dat");
  ParticleTransmissionFunc::SetAlphaFormalism(0);
  ParticleTransmissionFunc::SetNeutronFormalism(0);
  ParticleTransmissionFunc::SetProtonFormalism(0);
  GammaTransmissionFunc::SetEGDRType(1);
  GammaTransmissionFunc::SetPorterThomas(false);
  ParticleTransmissionFunc::SetPorterThomas(false);
  gsl_set_error_handler (&CoulFunc::GSLErrorHandler);
}
