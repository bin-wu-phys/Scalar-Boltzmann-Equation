/*
  Modified on 2015-04-16;
 */
#include <sstream>
#include "Init.h"
using namespace std;

class DSMC{
 private:
  double _msq;//mass squared of the particles
  int _n;//numbers of simulated particles (n) and of particles used for f (nf)
  //int _npair;//n(n-1)/2
  double _rho;//number density

  Init *_ran;//random number generator for momentum

  double **_pVec;//momentum of simulated particles

  double _Yhat;//max value of Y
  double _fmax;//the maximum value of f  

  double _lambda;//effective coupling: g^4/(64*pi) 

  double _dt;//time step
  double _compTime;//computation time

  ostringstream _fName;

  //dealing with the distribution function f
  int _np;//number of grids to contruct f
  double _pmin;//min momentum in the construction of f
  double _pmax;//max momentum in the construction of f

  double *_f,//constructed distribution
    *_p,//momentum for f
    *_vol,//volume element in momentum space to construct f
    _dp;//momentum increment

  double _nc;//number density of the condensate
  double _Vc;//Volume for condensed test particles 
  int _Nc;//number of condensated simulated particles

  //double _tRN;//time for generating random numbers
  //double _tScat;//time for computing scattering

public:
  double _t;
  DSMC(int n, int np,double pmin, double pmax, double rho, double t, double mass);
  ~DSMC();

  void getf();//construct f from simulated particles
  double getfMax();//get the maximal value of f
  void outputf();//output f;
  void output();//output

  double getYmax(int i, int j, bool);//give max Y given particles i and j;
  double getY(int i, int j, double *Omega);//get Y given particles i, j and Phi;
  void getYhat();//get Yhat
  double getE(int i);//get the energy of ith particle
  void nextTime();//time march
  void outputMac();//output macroscopic quantities
  void collision();//choose one pair to scatter
  void getPair(int &,int&);//get a pair of simulated particles
  void getOutgoingPair(int,int, double *,double *, double *);//get the outgoing momenta
  void updatePair(int,int,double *, double *);//update the momenta of the pair of simulated particles
  void updateDt();//calculate _dt
  double gettNext();
  void compTime();//calculate computation time
  void getnc();//get number density of the condensate
  void getNc();//get number of the condensated test particles

  void updatefin(double *p);//update the change of f due to the ingoing particles
  void updatefout(double *p);//update the change of f due to the outgoing particles
};
