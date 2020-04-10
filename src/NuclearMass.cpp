#include "NuclearMass.h"
#include "Constants.h"
#include <sstream>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <stdexcept>

#define HAS_EXP_MASS 1
#define HAS_TH_MASS  2

void NuclearMass::InitializeElements() {
#define ELEMENT(Z,EL) {elementTable_[std::string(EL)]=Z;}
#include "elements.h"
}

/**
 * @brief
 */
void NuclearMass::InitializeMasses(std::string filename) {
#ifndef MPI_BUILD
  std::cout << "Reading atomic masses..." << std::endl;
#endif
  std::ifstream in(filename.c_str());
  if(!in) {
    std::cout << "Could not read atomic masses data file." << std::endl;
    exit(1);
  }
  std::string line;
  for(int i = 0;i<4;++i) std::getline(in,line);
  while(!in.eof()) {
    std::getline(in,line);
    if(in.eof()) continue;
    std::istringstream lineStream(line);
    unsigned int mask = 0;
    int Z,A; lineStream>>Z>>A;
//  double atomicBindingE=(14.4381*pow(double(Z),2.39)+1.55468e-6*pow(double(Z),5.35))/1.e6;
    lineStream.ignore(5);
    char expMassExcessString[11];
    lineStream.read(expMassExcessString,10);
    expMassExcessString[10]='\0';
    std::istringstream expMassExcessStream(expMassExcessString);
    double expMassExcess;
    if(!(expMassExcessStream>>expMassExcess)) expMassExcess=0.;
    else mask |= HAS_EXP_MASS;
//  double expMass = (expMassExcess!=0.) ? A*uconv+expMassExcess-Z*eMass*uconv/1.e6+atomicBindingE : 0.;
    double expMass = (!!(mask&HAS_EXP_MASS)) ? A*uconv+expMassExcess : 0.;
    lineStream.ignore(10);
    char thMassExcessString[11];
    lineStream.read(thMassExcessString,10);
    thMassExcessString[10]='\0';
    std::istringstream thMassExcessStream(thMassExcessString);
    double thMassExcess; 
    if(!(thMassExcessStream>>thMassExcess)) thMassExcess=0.;
    else mask |= HAS_TH_MASS;
//  double thMass = (thMassExcess !=0.) ? A*uconv+thMassExcess-Z*eMass*uconv/1.e6+atomicBindingE : 0.;
    double thMass = (!!(mask&HAS_TH_MASS)) ? A*uconv+thMassExcess : 0.;
    char microCorrString[11];
    lineStream.read(microCorrString,10);
    microCorrString[10]='\0';
    std::istringstream microCorrStream(microCorrString);
    double microCorr; 
    if(!(microCorrStream>>microCorr)) microCorr=0.;
    massTable_[MassKey(Z,A)]=MassEntry(expMass,thMass,microCorr,mask);
  }
  in.close();
}

int NuclearMass::FindZ(std::string element) {
  ElementTable::const_iterator it = elementTable_.find(element);
  return (it!=elementTable_.end()) ? it->second : -1;
}

std::string NuclearMass::FindElement(int Z) {
  std::string returnString("?");
  for(ElementTable::const_iterator it = elementTable_.begin();
      it!=elementTable_.end();it++) {
  	if(it->second==Z) {
  	  returnString = it->first;
  	  break;
  	}
  }
  return returnString;
}

bool NuclearMass::FindMass(int Z, int A, double& M) {
  M = 0.;
  MassTable::const_iterator it = massTable_.find(MassKey(Z,A));
  if(it!=massTable_.end()) {
    if(!!(it->second.mask_&HAS_EXP_MASS)) M=it->second.expMass_;
    else if(!!(it->second.mask_&HAS_TH_MASS)) M=it->second.thMass_;
    else return false;
    return true;
  } else return false;
}


bool NuclearMass::MassDifference(int Z1, int A1, 
				 int Z2, int A2,
				 double& difference) {
  difference = 0.;
  double M1,M2;
  if(!FindMass(Z1,A1,M1)) return false;
  if(!FindMass(Z2,A2,M2)) return false;
  difference = M1-M2;
  return true;
}

double NuclearMass::MassDifference(int Z1, int A1, 
				 int Z2, int A2) {
  double M1,M2;
  if(!FindMass(Z1,A1,M1)) throw std::invalid_argument( "Cannot find M1");
  if(!FindMass(Z2,A2,M2)) throw std::invalid_argument( "Cannot find M1");
    return difference = M1-M2;
}

bool NuclearMass::QValue(int Z1, int A1,
			 int Z2, int A2,
			 double &qValue) {
  double massDifference;
  double decayMass;
  if(!MassDifference(Z1,A1,Z2,A2,massDifference)) return false;
  if(!FindMass(Z1-Z2,A1-A2,decayMass)) return false;
  qValue = massDifference-decayMass;
  return true;
}

bool NuclearMass::NeutronPairingGap(int Z, int A, double& pairingGap) {
  double massDifference1,massDifference2;
  if(!MassDifference(Z,A+1,Z,A,massDifference1)) return false;
  if(!MassDifference(Z,A-1,Z,A,massDifference2)) return false;
  pairingGap = (massDifference1+massDifference2)/2.;
  return true;
}

bool NuclearMass::ProtonPairingGap(int Z, int A, double& pairingGap) {
  double massDifference1,massDifference2;
  if(!MassDifference(Z+1,A+1,Z,A,massDifference1)) return false;
  if(!MassDifference(Z-1,A-1,Z,A,massDifference2)) return false;
  pairingGap = (massDifference1+massDifference2)/2.;
  return true;
}

bool NuclearMass::MicroEnergyCorr(int Z, int A, double& correction) {
  correction = 0.;
  MassTable::const_iterator it = massTable_.find(MassKey(Z,A));
  if(it!=massTable_.end()) {
    if(it->second.microEnergyCorr_!=0.) {
      correction = it->second.microEnergyCorr_;
      if(it->second.expMass_!=0.) {
	correction -= it->second.thMass_;
	correction += it->second.expMass_;
      }
      return true;
    } else return false;
  } else return false;
}

bool NuclearMass::HighestBoundEnergy(int Z, int A, double& energy) {
  double qValue[3];
  if(!QValue(Z,A,Z,A-1,qValue[0])) return false;
  if(!QValue(Z,A,Z-1,A-1,qValue[1])) return false;
  if(!QValue(Z,A,Z-2,A-4,qValue[2])) return false;
  energy=1000.;
  for(int i = 0;i<3;i++) 
    if(-1.*qValue[i]<energy) energy=-1.*qValue[i];
  return (energy!=1000.) ? true : false;
}

/*!
 *  Calculates liquid drop model mass based on TALYS parametrization.
 */

double NuclearMass::CalculateLDMMass(int Z, int A) {
	int N = A-Z;
	double dZ = double(Z);
	double dA = double(A);
	double dN = double(A-Z);
	double Mn = 8.07144;
	double Mp = 7.28899;
	double c1 = 15.677*(1.-1.79*pow((dN-dZ)/dA,2.));
	double c2 = 18.56*(1.-1.79*pow((dN-dZ)/dA,2.));
	double c3 = 0.717;
	double c4 = 1.21129;
	double del = 0.;
	if(Z%2==0&&N%2==0) del = -11./sqrt(dA);
	else if(Z%2!=0&&N%2!=0) del = 11./sqrt(dA);
	double mass = Mn*dN+Mp*dZ-c1*dA+c2*pow(dA,0.66666667)+c3*dZ*dZ/pow(dA,0.33333333)-
		c4*dZ*dZ/dA+del;
	return dA*uconv+mass;
}
