/**
 * @file Databases/Databases/NuclearMass.h
 * @author Mary Beard, Philipp Scholz
 * @date 2020
 * @brief It contains classes and methods for handling the nuclear masses: MassKey, MassEntry, MassTable, ElementTable, and NuclearMass
 */

#pragma once
#include <tr1/unordered_map>
#include <string>

/**
 * @brief Class for a key in a mass table
 */
class MassKey {
 public:
  MassKey(int Z, int A) : Z_(Z),A_(A) {};         /**< Constructor */
  bool operator<(const MassKey &right) const {
    if ( Z_ == right.Z_ ) return A_ < right.A_;
    else return Z_< right.Z_;                     /**< Definition of an operator for MassKey */
  };
  int Z_;                                         /**< Nuclear charge number */
  int A_;                                         /**< Nuclear mass number */
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

/**
 * @brief Class for an entry in a mass table.
 */
class MassEntry {
 public:
  /**
  * @brief Constructor of a MassEntry. Sets everything to 0.
  */
  MassEntry() {expMass_=0.;thMass_=0.;microEnergyCorr_=0.;mask_=0;};

  /**
   * @brief Constructor of a MassEntry, which initializes the member variables with the parameters.
   * @param expMass Experimental or recommended atomic mass excess in MeV of G. Audi, A.H. Wapstra, and C. Thibault (2003) Nucl. Phys. A **729**, 337
   * @param thMass Theoretical mass calculated FRDM atomic mass excess in MeV
   * @param microEnergyCorr The microscopic correction Emic corresponds to the difference between the total 
binding energy and the spherical macroscopic (droplet) energy.
   * @param mask Defines if this entry has a exp or th mass. 1 = has exp mass; 2 = has theo mass
   */
  MassEntry(double expMass,double thMass, double microEnergyCorr, unsigned int mask) :
  expMass_(expMass), thMass_(thMass), microEnergyCorr_(microEnergyCorr), mask_(mask) {};

  double expMass_; /**<Masses from experimental or recommended atomic mass excess in MeV of G. Audi, A.H. Wapstra, and C. Thibault (2003) Nucl. Phys. A **729**, 337*/
  double thMass_; /**< Masses from calculated FRDM atomic mass excess in MeV*/
  double microEnergyCorr_; /**< The microscopic correction Emic corresponds to the difference between the total binding energy and the spherical macroscopic (droplet) energy.*/
  unsigned int mask_; /**< Defines if this entry has a exp or th mass. 1 = has exp mass; 2 = has theo mass*/
};

typedef std::tr1::unordered_map<MassKey, MassEntry> MassTable; 
typedef std::tr1::unordered_map<std::string, int > ElementTable;

/**
 * @brief Class for handling masses of isotopes
 */
class NuclearMass {
 public:

  /**
   * @brief Initializes the map elementTable_
   * @details
   * Using a preprocessor macro the definition of entries in the map elementTable are carried
   * out by including elements.h.
   */
  static void InitializeElements();

  /** 
   * @brief Initializes the map massTable_
   * @param filename Path to the masstable file
   * @details
   * 1. Open the *masses.dat* file from the *tables/* directory.
   *  - Currently the FRDM95 masses table from RIPL-3 is used.
   *  - See [mass-frdm95.readme](https://www-nds.iaea.org/RIPL-3/masses/mass-frdm95.readme) for details.
   * 2. Skip the first four lines
   * 3. Read in mass excesses and calculate masses -> create a entry in the massTable_
   */
  static void InitializeMasses(std::string filename);

  /**
   * @brief Find the charge number of an element by its symbol
   * @param element The symbol of an element.
   * @return Atomic number as integer.
   */
  static int  FindZ(std::string element);

  /**
   * @brief Find an element symbol by its atomic number.
   * @param Z Atomic number as integer.
   * @return The symbol of an element as std::string.
   */  
  static std::string FindElement(int Z);

  /**
   * @brief Find the mass of an isotope in the massTable.
   * @param Z Atomic/charge number as integer.
   * @param A Mass number as integer.
   * @param M The mass as reference to a double.
   * @return True or false depending on whether the mass could be found or not.
   */
  static bool FindMass(int Z, int A, double& M);

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

  /**
   * @brief Calculate the QValue via the mass difference (M1-M2) and the mass of the decay product
   * @param Z1 Atomic number of isotope 1 (initial) as integer
   * @param A1 Mass number of isotope 1 (initial) as integer
   * @param Z2 Atomic number of isotope 2 (final) as integer
   * @param A2 Mass number of isotope 2 (final) as integer
   * @param qValue QValue for the decay as reference to a double
   * @return True or false depending on whether the qValue could be calculated or not.
   */
  static bool QValue(int Z1,int A1,int Z2,int A2,double& qValue);

  /**
   * @brief Calculates the Neutron Pairing GAp from the mass differences.
   * @param Z Atomic number as integer
   * @param A Mass number as integer
   * @param pairingGap Pairing gap as reference to a double.
   * @return True or false depending on whether the pairingGap could be calculated
   */
  static bool NeutronPairingGap(int Z,int A,double& pairingGap);

  /**
   * @brief Calculates the Proton Pairing GAp from the mass differences.
   * @param Z Atomic number as integer
   * @param A Mass number as integer
   * @param pairingGap Pairing gap as reference to a double.
   * @return True or false depending on whether the pairingGap could be calculated
   */
  static bool ProtonPairingGap(int Z,int A,double& pairingGap);

  /**
   * @brief Finds the micro energy correction from the massTable_.
   * @param Z Atomic number as integer.
   * @param A Mass number as integer.
   * @param correction Reference to a double for the micro energy correction.
   * @return True or false depending on whether the micro energy correction could be found or not.
   */
  static bool MicroEnergyCorr(int Z,int A,double& correction);

  /**
   * @brief Calculates the qValues for neutron, proton and alpha decay and puts the negative qValue as energy.
   * @param Z Atomic number as integer.
   * @param A Mass number as integer.
   * @param energy The energy of the highest bound state as negative qValue.
   * @return True or false whether an energy could be set or not.
   */
  static bool HighestBoundEnergy(int Z,int A, double& energy);

  /**
   * @brief Calculates liquid drop model mass based on TALYS parametrization.
   * @param Z Atomic number as integer.
   * @param A Mass number as integer.
   * @return The LDM mass as double.
   */
  static double CalculateLDMMass(int Z,int A);

 private:
  static MassTable massTable_; /**< Map of MassKey and MassEntry*/
  static ElementTable elementTable_; /**< Map of atomic number and element symbol*/
};

