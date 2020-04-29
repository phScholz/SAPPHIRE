/**
 * @file LevelDensityTable.h
 * @date 2020-04-27
 * @brief Declaration and Definition of the LevelDensityTable class
 */

#pragma once
#include "LevelDensity.h"
#include <string>
#include <vector>
#include <iostream>
#include <tr1/unordered_map>

class LevelDensityTable : public LevelDensity{
    
    protected:
        virtual void GetFileName() = 0;
        virtual void ReadFile() = 0;
        virtual void FillVector() = 0;
        virtual double TotalLevelDensity(double E) = 0;
        virtual double CalculateDensity(double E) = 0;
        void SetFilename(std::string x){filename=x;}

    public:
        LevelDensityTable(int Z, int A, double J, int parity): LevelDensity(Z,A,J,parity){

        }

    protected:
        std::vector< std::pair<double, double> > DensityVector; /**< Densities will be read to this vector*/
        std::string filename; /**< Filename of the data table*/
};

class ChargeMassParityKey {
 public:
  ChargeMassParityKey(int Z, int A, int P) : Z_(Z),A_(A), P_(P){};    /**< Constructor */

  bool operator <(const ChargeMassParityKey &right) const {
    if ( Z_ == right.Z_ ) return A_ < right.A_;
    else return Z_< right.Z_;                     /**< Definition of an operator for MassKey */
    if (Z_== right.Z_ && A_==right.A_) return P_ < right.P_;
  };

  int Z_;                                         /**< Nuclear charge number */
  int A_;                                         /**< Nuclear mass number */
  int P_;                                          /**< Parity +1 or -1*/
};

namespace std {
  namespace tr1 {
    template<>
    struct hash<ChargeMassParityKey> {
      std::size_t operator() (ChargeMassParityKey const &key) const {
	size_t hash = 23;
	hash = (hash*37) + key.Z_;
	hash = (hash*37) + key.A_;
    hash = (hash*37) + key.P_;
	return hash;
      }
    };
  }

  template<>
  struct equal_to<ChargeMassParityKey> {
    bool operator()(ChargeMassParityKey const &left, ChargeMassParityKey const &right) const {
      if ( left.Z_ == right.Z_ && left.P_ == right.P_) return left.A_ == right.A_;
      else return false;
    }
  };
}
