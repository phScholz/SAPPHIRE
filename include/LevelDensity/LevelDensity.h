/**
 * @file LevelDensity.h
 * @date 2020-04-29
 * @brief Declaration and definition of LevelDensity class and it's members
 */
#pragma once
#include <string>

/**
 * @brief Base class for all LevelDensity models.
 * @details
 * Important are 
 * - TotalLevelDensity()
 * - CalculateDensity()
 * 
 * These methods should be provided by all other NLD models.
 */
class LevelDensity {
 public:
 /**
  * @brief Constuctor of a LevelDensity.
  * @param Z Charge number
  * @param A Mass number
  * @param J Spin
  * @param P Parity
  */
 LevelDensity(int Z, int A, double J, int P) :
  Z_(Z), A_(A), J_(J), parity_(P) {
  };

  /**
   * @brief Destructor
   */
  virtual ~LevelDensity(){};

  /**
   * @brief Operator which is used sometimes in the code to call CalculateDensity().
   */
  double operator()(double E){
    return CalculateDensity(E);
  };
  
  /**
   * @brief For some reasons, these classes need to have access to some functions.
   */
  friend class KopeckyUhlGSF;
  friend class TransitionRateFunc;

 protected:
  /**
   * @brief Method to set tables
   */ 
  void SetTables(bool x){tables = x;};

  /**
   * @brief Method to calculate the total Level Density.
   */
  virtual double TotalLevelDensity(double E){};

  /**
   * @brief Method to calculate the partial density.
   * @param E Excitation Energy. 
   */
  virtual double CalculateDensity(double E){};

  

 protected:

  int Z_; /**< Charge number of the nucleus*/
  int A_; /**< Mass number of the nucleus*/
  int parity_; /**< Parity of the states*/
  double J_; /**< Spin of the states*/
  bool tables = false; /**< Boolean to decide if its a LevelDensityTable or Formular.*/
};


