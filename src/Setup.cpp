#include "NuclearMass.h"
#include "GammaTransmissionFunc.h"
#include "NuclearLevels.h"
#include "Decayer.h"
#include "Sapphire_config.h"
#include "TransitionRateFunc.h"
#ifndef MPI_BUILD
#include "CrossSection.h"
#endif
#include "PreEqDecayer.h"
#include "ParticleTransmissionFunc.h"
#include "CoulFunc.h"
#include <iostream>
#include <gsl/gsl_errno.h>

#ifndef MPI_BUILD
bool CrossSection::residualGamma_;
bool CrossSection::residualNeutron_;
bool CrossSection::residualProton_;
bool CrossSection::residualAlpha_;
bool CrossSection::calculateGammaCutoff_;
std::vector<double> CrossSection::rateTemps_;
std::vector<double> CrossSection::macsEnergies_;
#endif
bool Decayer::isCrossSection_;
bool PreEqDecayer::isCrossSection_;
double Decayer::maxL_;
double PreEqDecayer::maxL_;
double TransitionRateFunc::gammaCutoffEnergy_;
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
#ifndef MPI_BUILD
  CrossSection::SetResidualGamma(true);
  CrossSection::SetResidualNeutron(false);
  CrossSection::SetResidualProton(false);
  CrossSection::SetResidualAlpha(false);
  CrossSection::SetCalculateGammaCutoff(true);
  CrossSection::CreateTempVector();
  CrossSection::CreateMACSEnergiesVector();
#endif
  Decayer::SetCrossSection(false);
  PreEqDecayer::SetCrossSection(false);
  Decayer::SetMaxL(8.);
  PreEqDecayer::SetMaxL(8.);
  TransitionRateFunc::SetGammaCutoffEnergy(10000.);
  NuclearMass::InitializeElements();
  NuclearMass::InitializeMasses(sourceDirectory()+"/tables/masses.dat");
  GammaTransmissionFunc::InitializeGDRParameters(sourceDirectory()+"/tables/ripl3_gdr_parameters.dat");
  NuclearLevels::InitializeLevels(sourceDirectory()+"/levels/",sourceDirectory()+"/tables/spinod.dat");
  ParticleTransmissionFunc::SetAlphaFormalism(0);
  ParticleTransmissionFunc::SetNeutronFormalism(0);
  ParticleTransmissionFunc::SetProtonFormalism(0);
  GammaTransmissionFunc::SetEGDRType(1);
  GammaTransmissionFunc::SetPorterThomas(false);
  ParticleTransmissionFunc::SetPorterThomas(false);
  gsl_set_error_handler (&CoulFunc::GSLErrorHandler);
}
