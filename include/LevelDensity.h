#pragma once
#include <string>
class LevelDensity {
 public:
 LevelDensity(int Z, int A, double J, int parity) :
  Z_(Z), A_(A), J_(J), parity_(parity) {
  };
  virtual ~LevelDensity(){};
  double operator()(double E){
    return CalculateDensity(E);
  };
  
  friend class KopeckyUhlGSF;
  friend class TransitionRateFunc;

 protected:
  
  
  void SetTables(bool x){tables = x;};

  virtual double TotalLevelDensity(double E){};
  virtual double CalculateDensity(double E){};

  

 protected:

  int Z_; /**< Charge number of the nucleus*/
  int A_; /**< Mass number of the nucleus*/
  int parity_; /**< Parity of the states*/
  double J_; /**< Spin of the states*/
  bool tables = false;
};


