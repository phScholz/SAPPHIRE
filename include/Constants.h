#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <complex>
#include <vector>
#include <cstdlib>

const double pi=3.141592650;
const double hbarc=197.32696310;
const double uconv=931.4940880;
const double fstruc=1.00/137.0359996790;
const double boltzConst=8.6171e-2;
const double lightSpeedInCmPerS=29979245800.;
const double avagadroNum=6.02214179e23;
const double eMass = 548.579894;

typedef std::complex<double> complex;
typedef std::vector<double> vector_r;
typedef std::vector<std::complex<double> > vector_c;
typedef std::vector<std::vector<double> > matrix_r;
typedef std::vector<std::vector<std::complex<double> > > matrix_c;
typedef std::vector<std::vector<std::vector<double> > > vector_matrix_r;
typedef std::vector<std::vector<std::vector<std::complex<double> > > > vector_matrix_c;

#endif