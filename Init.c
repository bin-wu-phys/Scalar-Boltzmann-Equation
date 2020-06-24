/*
  class Ranx generates the random number for the longitudianl fraction of the energy carried by the child gluons. Created on 2014-09-20 by Bin Wu.
  Modified on 2015-04-09, add |p| to generate.
 */

#include "Init.h"
#include <math.h>
#include <iostream>
#include "const.h"

//using namespace std;



Init::Init(Ullong seeds, double cmax){
  _seeds=seeds;_cmax=cmax;
  _ran=new Ran(_seeds);
}

Init::~Init(){
  delete _ran;
}

double Init::doub(){
  return _ran->doub();
}

void Init::generate(double *x){//f_0=3*theta(1-p)/(4*pi*cmax^3)
  double c3=_ran->doub(), c=_cmax*pow(c3,1.0/3.0), ct=2.0*(_ran->doub()-0.5),st=sqrt(1.0-ct*ct),phi=2.0*PI*_ran->doub();
  x[0]=c*st*cos(phi);x[1]=c*st*sin(phi);x[2]=c*ct;x[3]=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
}

void Init::cgc(double *x){//f_0=theta(cmax-p_T)/(4*pi*cmax^2)
  double c2=_ran->doub(), c=_cmax*sqrt(c2), phi=2.0*PI*_ran->doub();
  x[0]=c*cos(phi);x[1]=c*sin(phi);x[2]=0.0;
}

void Init::sphere(double *x){
  double ct=2.0*(_ran->doub()-0.5),st=sqrt(1.0-ct*ct),phi=2*PI*_ran->doub();
  x[0]=st*cos(phi);x[1]=st*sin(phi);x[2]=ct;
}

int Init::ranInt(int N){
  int rInt;
  rInt=floor(_ran->doub()*N);
  if(rInt==N) rInt-=1;
  return rInt;
}

void Init::save(ofstream *out){//save parameters
  out->write(reinterpret_cast<char*>(&_ran->u),(sizeof _ran->u));
  out->write(reinterpret_cast<char*>(&_ran->v),(sizeof _ran->v));
  out->write(reinterpret_cast<char*>(&_ran->w),(sizeof _ran->w));
  out->write(reinterpret_cast<char*>(&_cmax),(sizeof _cmax));
  out->write(reinterpret_cast<char*>(&_seeds),(sizeof _seeds));
}

void Init::read(ifstream *in){//read parameters
  in->read(reinterpret_cast<char*>(&_ran->u),(sizeof _ran->u));
  in->read(reinterpret_cast<char*>(&_ran->v),(sizeof _ran->v));
  in->read(reinterpret_cast<char*>(&_ran->w),(sizeof _ran->w));
  in->read(reinterpret_cast<char*>(&_cmax),(sizeof _cmax));
  in->read(reinterpret_cast<char*>(&_seeds),(sizeof _seeds));
}
