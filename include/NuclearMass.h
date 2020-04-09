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
  static bool MassDifference(int,int,int,int,double&);
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
