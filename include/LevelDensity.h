#pragma once
#include <string>
class LevelDensity {
 public:
 LevelDensity(int Z, int A, double J) :
  Z_(Z), A_(A), J_(J) {
    criticalU_=2.5+150./A;
  };
  virtual ~LevelDensity() {};
  double operator()(double);
  double TotalLevelDensity(double);
  double CalculateDensity(double E);
  double LookUpDensity(double E);
  friend class KopeckyUhlGSF;

 protected:
  virtual void CalcBackShift() = 0;
  virtual double CalcDensityParam(double) = 0;
  virtual double CalcNuclearTemp(double) = 0;
  //virtual void ReadLevelDensityFile(std::string fileName);
  //virtual double InterpolateDensity(double energy);
  void SetTables(bool x){tables = x;};
  //virtual std::string GetLevelDensityFile();
  void CalcConstantTempTerms();

 protected:
  static constexpr double zeta_ = 1.0;
  static constexpr double r0_ = 1.25;
  int Z_;
  int A_;
  double J_;
  double backshift_;
  double criticalU_;
  double constAngTerm_;
  double nuclearTemp_;
  double e0_;
  bool tables = false;
  std::string file;
};


