#include "EquivSquareWell.h"

double EquivSquareWell::CalcTransmission(double s, int l, double energy) {
  double rho = sqrt(2.*uconv)/hbarc*GetR1()*sqrt(redmass_*energy);
  double pene = coulFunc_->Penetrability(l,GetR1(),energy);
  return 1-exp(-4.*pi*pene*GetF()/pi/rho/sqrt(1+GetV1()/energy));
}

double EquivSquareWell::GetF() {
  double f=0.;
  if(pType_==0) {
    f = 2.7;
  } else if(pType_==1) {
    f = 2.7;
  } else if(pType_==2) {
    f = 4.8;
  } 
  return f;
}

double EquivSquareWell::GetV0() {
  double v=0.;
  if(pType_==0) {
    v = 50;
  } else if(pType_==1) {
    v = 50;
  } else if(pType_==2) {
    v = 60;
  } 
  return v;
}

double EquivSquareWell::GetV1() {
  return GetV0()*pow(GetR0()/GetR1(),2.);
}

double EquivSquareWell::GetR0() {
  double radius=0.;
  if(pType_==0) {
    radius = 1.25*pow(double(m2_),0.3333333333);
  } else if(pType_==1) {
    radius = 1.25*pow(double(m2_),0.3333333333);
  } else if(pType_==2) {
    radius = 1.09*pow(double(m2_),0.3333333333)+1.6;
  } 
  return radius;
}

double EquivSquareWell::GetR1() {
  double radius=0.;
  if(pType_==0) {
    radius = GetR0()+0.1;
  } else if(pType_==1) {
    radius = GetR0()+0.1;
  } else if(pType_==2) {
    radius = GetR0()+0.7;
  } 
  return radius;
}
