#ifndef MCFADDENSATCHLERPOTENTIAL_H
#define MCFADDENSATCHLERPOTENTIAL_H

#include "Potential.h"

class TransmissionFunc;

class McFaddenSatchlerPotential : public Potential {
 public:
  McFaddenSatchlerPotential(int,int,int,int,double,int,
			    double,int,double,int,double,
			    double,double,double,
			    TransmissionFunc*); 

  std::complex<double> Calculate(double,int,double,double,double) const;
 private:
    double rR_;
    double rI_;
    double aR_;
    double aI_;
    double V_;
    double W_;
};

#endif
