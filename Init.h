/*
  class Ranx generates the random number for the longitudianl fraction of the energy carried by the child gluons. Created on 2014-09-20 by Bin Wu.
 */
#include <fstream>
#include "ran.h"
using namespace std;

class Init{
 private:
  Ran *_ran;
  double _cmax;//max velocity
  Ullong _seeds;

 public:
  Init(Ullong seeds, double cmax);
  ~Init();
  
  double doub();
  void generate(double *x);//f_0=3*theta(1-p)/(4*pi*cmax^3)
  void cgc(double *x);//f_0=theta(cmax-p_T)/(4*pi*cmax^2)
  void sphere(double *x);
  int ranInt(int N);
  void save(ofstream*);
  void read(ifstream*);
};
