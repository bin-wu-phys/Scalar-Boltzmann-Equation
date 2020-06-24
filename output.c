/*
  output methods of class DSMC
 */

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "dsmc.h"
#include <cmath>
#include <omp.h>
using namespace std;

double DSMC::gettNext(){//get next time to save
  int ex=floor(log10(_t));ex--;
  int rt=ceil(_t*pow(10.0,-(double)ex));
  //cout << _t << ": " << ex << " " << rt << endl; 
  return ((double)rt)*pow(10.0,(double) ex);
}

void DSMC::outputf(){
  ofstream out;
  out.precision(15);
  out.open(string(_fName.str()+".dat").c_str(),ios::app);

  out << "# t = " << _t << endl;
  
  for(int i=0;i<_np;i++){
    out << _p[i] << " " << _f[i] << endl;
  }
  out <<"\n\n";
  out.close();
}

void DSMC::outputMac(){//output macroscopic quantities

  double PT=0,PL=0,eps=0.0,p4=0.0, Ep, pz2, pT2, p2=0.0;
  ofstream out;
  out.precision(15);
  out.open(string(_fName.str()+"Mac.dat").c_str(),ios::app);
  getnc();
  for(int i=0;i<_n;i++){
    pT2=_pVec[i][0]*_pVec[i][0]+_pVec[i][1]*_pVec[i][1];
    pz2=_pVec[i][2]*_pVec[i][2];p2+=pT2+pz2;p4+=((pT2+pz2)*(pT2+pz2));
    Ep=sqrt(pT2+pz2+_msq);
    PL+=(pz2/Ep);PT+=(pT2/Ep);eps+=Ep;
  }
  out << _t << " " << _rho*eps/((double)_n) << " " <<  0.5*_rho*PT/((double)_n) << " " <<  _rho*PL/((double)_n)  << " " << _rho*p2/((double)_n) << " " << _rho*p4/((double)_n) << " " << _nc << " " << _Nc << " " << _dt << " " << omp_get_wtime() - _compTime << endl;
  out.close();
}

void DSMC::output(){
  outputMac();outputf();
}

void printTime(ofstream &out, double _compTime){
  if(_compTime>=60.0){
    if(_compTime>=3600){
       out << _compTime/3600.0 << " hours.\n\n\n";
    }
    else {
      out <<  _compTime/60.0 << " minutes.\n\n\n";
    }
  }
  else {
    out << _compTime << " seconds.\n\n\n";
  }
}

void DSMC::compTime(){
  //Compuation time
  _compTime= omp_get_wtime() - _compTime;
  ofstream out;out.precision(15);
  out.open(string(_fName.str()+"Mac.dat").c_str(),ios::app);
  out << "#Calculation is done and the running time is "; printTime(out,_compTime);
  //out << "#Time for generateing random numbers is "; printTime(out,_tRN);
  //out << "#Time for computing scattering is "; printTime(out,_tScat);
  out.close();
}

