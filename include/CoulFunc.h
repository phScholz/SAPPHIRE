#ifndef COULFUNC_H
#define COULFUNC_H
#pragma once
struct CoulWaves {
  double F; 
  double dF;
  double G;
  double dG;
};

class CoulFunc {
 public:
  CoulFunc(int z1, int z2, double redmass, bool useGSLFunctions);
  int z1() const;
  int z2() const;
  double redmass() const;
  int lLast() const;
  double radiusLast() const;
  double energyLast() const;
  struct CoulWaves coulLast() const;
  void setLast(int, double, double, CoulWaves);
  CoulWaves operator()(int,double,double);
  double Penetrability(int,double,double);
  double PEShift(int,double,double);
  double PEShift_dE(int,double,double);;
  static void GSLErrorHandler(const char*, const char*, int, int);
 private:
  static double thisPEShift(double,void*);
  typedef struct DEShiftParams {
    CoulFunc *coulFunc;
    int lValue;
    double radius;
  } DEShiftParams;
  DEShiftParams dEShiftParams_;
  bool useGSLFunctions_;
  int z1_;
  int z2_;
  int lLast_;
  double redmass_;
  double radiusLast_;
  double energyLast_;
  struct CoulWaves coulLast_;
};

#endif
