/**
 * @file NuclearMass.h
 * 
 */

#ifndef NUCLEARMASS_H
#define NUCLEARMASS_H

#include <tr1/unordered_map>
#include <string>

class MassKey {
 public:
  MassKey(int Z, int A) : Z_(Z),A_(A) {};
  bool operator<(const MassKey &right) const {
    if ( Z_ == right.Z_ ) return A_ < right.A_;
    else return Z_< right.Z_;
  };
  int Z_;
  int A_;
};

namespace std {
  namespace tr1 {
    template<>
    struct hash<MassKey> {
      std::size_t operator() (MassKey const &key) const {
	size_t hash = 23;
	hash = (hash*37) + key.Z_;
	hash = (hash*37) + key.A_;
	return hash;
      }
    };
  }
  template<>
  struct equal_to<MassKey> {
    bool operator()(MassKey const &left, MassKey const &right) const {
      if ( left.Z_ == right.Z_ ) return left.A_ == right.A_;
      else return false;
    }
  };
}

class MassEntry {
 public:
  MassEntry() {expMass_=0.;thMass_=0.;microEnergyCorr_=0.;mask_=0;};
 MassEntry(double expMass,double thMass, double microEnergyCorr, unsigned int mask) :
  expMass_(expMass), thMass_(thMass), microEnergyCorr_(microEnergyCorr), mask_(mask) {};
  double expMass_;
  double thMass_;
  double microEnergyCorr_;
  unsigned int mask_;
};

typedef std::tr1::unordered_map<MassKey, MassEntry> MassTable;
typedef std::tr1::unordered_map<std::string, int > ElementTable;

class NuclearMass {
 public:
  static void InitializeElements();
  static void InitializeMasses(std::string);
  static int  FindZ(std::string);
  static std::string FindElement(int);
  static bool FindMass(int, int, double&);

/**
* @brief Calculate the mass difference between two isotopes and writes it in "difference". Returns False if the mass cannot be found.
* @param Z1 Charge of Isotope 1
* @param A1 Massnumber of Isotope 1
* @param Z2 Charge of Isotope 2
* @param A2 Massnumber of Isotope 2
* @param difference reference to the mass difference double
* @return Returns a *True* or *False* depending on if the masses of both isotopes have been found by FindMass()
*/
  static bool MassDifference(int,int,int,int,double&);

/**
* @brief Calculate the mass difference between two isotopes and returns it as double. If the mass cannot be found, an exception is thrown.
* @param Z1 Charge of Isotope 1
* @param A1 Massnumber of Isotope 1
* @param Z2 Charge of Isotope 2
* @param A2 Massnumber of Isotope 2
* @return Returns a *True* or *False* depending on if the masses of both isotopes have been found by FindMass()
*/
  static double MassDifference(int,int,int,int);


  static bool QValue(int,int,int,int,double&);
  static bool NeutronPairingGap(int,int,double&);
  static bool ProtonPairingGap(int,int,double&);
  static bool MicroEnergyCorr(int,int,double&);
  static bool HighestBoundEnergy(int,int,double&);
  static double CalculateLDMMass(int,int);
 private:
  static MassTable massTable_;
  static ElementTable elementTable_;
};

#endif
