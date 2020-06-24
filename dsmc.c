/*
  Modified on 2015-05-24: _dt changed.
 */
#include "const.h"
#include "dsmc.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <omp.h>
using namespace std;

DSMC::DSMC(int n, int np, double pmin, double pmax, double rho, double t, double mass){
  _compTime=omp_get_wtime();
  _n=n;_np=np;_pmax=pmax;_pmin=pmin;//_npair=_n*(_n-1)/2;
  _rho=rho;_t=t;_lambda=1.0;
  _ran=new Init(11,1.0);_msq=mass*mass;
  //_tRN=0.0;_tScat=0.0;

  //cout << "f_0 = " << (6.0*PI*PI)*_rho  << endl;
  _pVec=new double *[_n];for(int i=0;i<_n;i++) _pVec[i]=new double [4];

  //get momenta of test particles
  for(int i=0;i<_n;i++){
    _ran->generate(_pVec[i]); if(_pVec[i][3]<_pmin) _Nc++;
  }
  //cout << "p_" << i << " = (" << _pVec[i][0] << ", "  << _pVec[i][1] << ", "  << _pVec[i][2] << ")" << endl;  
 
  _dp=(_pmax-_pmin)/((double)_np);

  _f=new double [_np];_p=new double [_np];_vol=new double [_np];

  double vol=_rho*(PI*PI)*6.0/((double)_n);
  _Vc=vol/(_pmin*_pmin*_pmin);
  //cout << "_Vc: " << _Vc <<endl;
  for(int i=0;i<_np;i++) {
    _p[i]=_pmin+(i+0.5)*_dp;
    _vol[i]=vol/((pow(_pmin+double(i+1)*_dp,3.0)-pow(_pmin+double(i)*_dp,3.0)));
  }

  //get the number of condensated simulated particles
  getNc();
  cout << "Nc:  " << _Nc << endl;
  //get f
  getf();_fmax=getfMax();

  //get Yhat
  getYhat();
  cout << "Yhat = " << _Yhat << endl;

  //get file name
  int nrho=floor(log10(_rho));
  if(nrho<0)
    _fName << "rho" << int(_rho*pow(10.0,-nrho))  << "em" << -nrho;
  else
     _fName << "rho" << int(_rho*pow(10.0,-nrho))  << "e" << nrho;
  _fName << "m"  << setfill('0') << setw(3) << int(100.0*mass);
  _fName << "N"  << log10(_n);
  _fName << "M"  << _np;
  int ip=-floor(log10(_pmin));
  _fName << "pmin" <<  int(_pmin*pow(10.0,ip)) <<"em"  << ip;

  _fName << "t0" << int(_t);
  _fName << "v3.1.4";

  //set time step
  updateDt();
  cout << "dt: " << _dt << endl;

  // current date/time based on current system
  time_t now = time(0);
  
  // convert now to string form
  char* date = ctime(&now);
  cout << "V3.1.4:\n" <<endl;
  cout << "The starting time is: " << date << endl;
  cout << "n = " << _n << endl; 
  cout << "m = " << _np << endl;
  cout << "rho = " << _rho << endl;
  cout << "pmin = " << _pmin << endl;
    
}



DSMC::~DSMC(){
  delete _ran;
  for(int i=0;i<_n;i++) delete [] _pVec[i]; delete [] _pVec;
  delete [] _f;delete [] _p;delete [] _vol;
}


void DSMC::updatefout(double *p){//update the change of f due to the outgoing scattering
  if(p[3]>=_pmin){
    int index=floor((p[3]-_pmin)/_dp);
    if(index<_np){
      _f[index]+=_vol[index];
      if(_f[index]>_fmax) _fmax=_f[index];
    }
  }
}

void DSMC::updatefin(double *p){//update the change of f due to the ingoing scattering
  if(p[3]>=_pmin){
    int index=floor((p[3]-_pmin)/_dp);
    if(index<_np){
      _f[index]-=_vol[index];
    }
  }
}


void DSMC::getf(){
  for(int i=0;i<_np;i++) _f[i]=0.0;
  for(int i=0;i<_n;i++){
    //updatef(_pVec[i]);
    if(_pVec[i][3]>=_pmin){
      int index=floor((_pVec[i][3]-_pmin)/_dp);
      if(index<_np) _f[index]+=_vol[index];
    }
  }
  /*  for(int i=0;i<_np;i++){
    _f[i]*=_vol[i];
    }*/
}

double DSMC::getfMax(){
  double fmax=0;
  for(int i=0;i<_np;i++){//discard the condensate
    if(_f[i]>fmax){fmax=_f[i];}
  }
  return fmax;
}


double DSMC::getE(int i){
  // cout << "E_p" << i << " = " << sqrt(psq+_msq) << endl;
  return sqrt(_pVec[i][3]*_pVec[i][3]+_msq);
}

double DSMC::getYmax(int i, int j, bool isFinal){
  double Ep=getE(i), Ep1=getE(j);
  double Etotsq=Ep+Ep1;Etotsq*=Etotsq;
  double Ptotsq=0;
  for(int d=0;d<3;d++){
    double pij=_pVec[i][d]+_pVec[j][d];
    Ptotsq+=pij*pij;
  }
  double s=Etotsq-Ptotsq;double fenh;//Bose enhancement
  fenh=1.0+_fmax;
  if(isFinal){
    fenh+=double(_Nc)*_Vc;
  }else{
    fenh+=_fmax;
  }

  return 1.2*Etotsq*sqrt(s-4.0*_msq)*fenh/(Ep*Ep1*pow(s,1.5));
}

void DSMC::getYhat(){//get Yhat, only at initial time without a condensate
    /*_Yhat=0;
   for(int i=0;i<_n;i++){
    for(int j=i+1;j<_n;j++){
    double Y=getYmax(i,j);
      if(Y>_Yhat) _Yhat=Y;
    }
  }
    */
  double fmax=1.0+2.0*_fmax;

  double msqInv=1.0/_msq;
  double emax=sqrt(1.0+msqInv);
  _Yhat=msqInv*(1.0+emax)*(1+emax)*fmax/(6.0*sqrt(3.0)*emax);
}

double DSMC::getY(int i, int j, double *Omega){//get Y given particles i, j, Phi, pp=pv, p1p=pw;
  double Ep=getE(i), Ep1=getE(j);
  //cout << "{ " << _pVec[i][0] << ", "  << _pVec[i][1] << ", " << _pVec[i][2] << "}" << endl;  
  //cout << "{ " << _pVec[j][0] << ", "  << _pVec[j][1] << ", " << _pVec[j][2] << "}" << endl;
  //cout << "{ " << Omega[0] << ", "  << Omega[1] << ", " << Omega[2] << "}" << endl;  

  double Etotsq=Ep+Ep1;Etotsq*=Etotsq;
  double Ptotsq=0,PtotOmega=0, Ptot[3];
  for(int d=0;d<3;d++){
    Ptot[d]=_pVec[i][d]+_pVec[j][d];
    PtotOmega+=Ptot[d]*Omega[d];
    Ptotsq+=Ptot[d]*Ptot[d];
  }
  double s=Etotsq-Ptotsq;
  return Etotsq*sqrt(s-4.0*_msq)/(Ep*Ep1*pow(Etotsq-PtotOmega*PtotOmega,1.5));
}


void DSMC::nextTime(){

  //_t+=_dt;
  collision();

}

void DSMC::collision(){//choose one pair to scatter

  //get a pair
  //double tRN=omp_get_wtime();
  int v,w; getPair(v,w);
  //double Omega[3];
  //_ran->sphere(Omega);
  //double x1=_ran->doub();
  //_tRN+=omp_get_wtime()-tRN;

  //tRN=omp_get_wtime();
  if(_pVec[v][3]>=_pmin||_pVec[w][3]>=_pmin){//at most one is a condensate
    //get Omega
    double Omega[3];
    _ran->sphere(Omega);
    
    double pv[4],pw[4];
    getOutgoingPair(v,w,Omega,pv,pw);
 
    if(pv[3]>=_pmin||pw[3]>=_pmin){//at most one is a condensate
      bool isInit=((_pVec[v][3]<_pmin||_pVec[w][3]<_pmin)&&pv[3]>=_pmin&&pw[3]>=_pmin);
      bool isFinal=((_pVec[v][3]>=_pmin&&_pVec[w][3]>=_pmin)&&(pv[3]<_pmin||pw[3]<_pmin));
      bool isPart=(_pVec[v][3]>=_pmin&&_pVec[w][3]>=_pmin&&pv[3]>=_pmin&&pw[3]>=_pmin);
      if(isInit||isFinal||isPart){
	_t+=_dt;
       //get Y
       double Y=1.0;
       
       int vp=floor((pv[3]-_pmin)/_dp);
       int vpbad=vp;
       if(vp<0) Y+=double(_Nc)*_Vc;
       else if(vp<_np) Y+=_f[vp];
       vp=floor((pw[3]-_pmin)/_dp);
       if(vp<0) Y+=double(_Nc)*_Vc;
       else if(vp<_np) Y+=_f[vp];
       Y= getY(v,w,Omega)*Y/_Yhat;

       if(Y>1.0){
	 cout << "(1+fmax)^2 = " << (1.0+getfMax())*(1.0+getfMax()) << endl;
	 cout << "f[" << vpbad << "] = " << _f[vpbad] << ", f[" << vp << "] = " << _f[vp] << endl; 
	 cout << "Yhat is not good = " << _Yhat << " and Y = " << Y << endl;
	 //exit(0);
       }  
       
       //If v, w need to be updated
       double x1=_ran->doub();
       if(Y>x1){
	 updatefin(_pVec[v]);updatefin(_pVec[w]);
	 updatefout(pv);updatefout(pw);updatePair(v,w,pv,pw);
	 if(isInit) _Nc--;
	 else if(isFinal) _Nc++;

	 double Ymax=getYmax(v,w,isFinal);
	 if(Ymax>_Yhat){
	   _Yhat=Ymax;updateDt();
	 }
       }
      }  
   }
  }
  //_tScat+=omp_get_wtime()-tRN;
}


void DSMC::getPair(int &v, int &w){
  do{
    v=_ran->ranInt(_n);w=_ran->ranInt(_n);
  }while(v>=w);
}

/*
void DSMC::getPair(int &v, int &w){
  int x =_ran->ranInt(_npair);
  v=floor(sqrt(2.0*x+2.25)-0.5);
  w=x+1-v*(v+1)/2;
  //cout << x << ": (" << v << ", " << w << ")\n";
}
*/

void DSMC::getOutgoingPair(int i, int j, double *Omega, double *pv, double *pw){
  double Ep=getE(i), Ep1=getE(j);
  
  double Etot=Ep+Ep1,Etotsq=Etot*Etot;
  double Ptotsq=0,PtotOmega=0, Ptot[3];
  for(int d=0;d<3;d++){
    Ptot[d]=_pVec[i][d]+_pVec[j][d];
    PtotOmega+=Ptot[d]*Omega[d];
    Ptotsq+=Ptot[d]*Ptot[d];
  }
  double s=Etotsq-Ptotsq;
  double a=Etot*sqrt((s-4.0*_msq)/(Etotsq-PtotOmega*PtotOmega));
  
  for(int d=0;d<3;d++){
    pv[d]=0.5*(Ptot[d]+a*Omega[d]);
    pw[d]=0.5*(Ptot[d]-a*Omega[d]);
  }
  pv[3]=sqrt(pv[0]*pv[0]+pv[1]*pv[1]+pv[2]*pv[2]);
  pw[3]=sqrt(pw[0]*pw[0]+pw[1]*pw[1]+pw[2]*pw[2]);
  //cout << "{ " << _pVec[i][0] << ", "  << _pVec[i][1] << ", " << _pVec[i][2] << "}" << endl;  
  //cout << "{ " << _pVec[j][0] << ", "  << _pVec[j][1] << ", " << _pVec[j][2] << "}" << endl;
}

void DSMC::updatePair(int i, int j, double *pv, double *pw){
  
  for(int d=0;d<4;d++){
    _pVec[i][d]=pv[d];
    _pVec[j][d]=pw[d];
  }
}


void DSMC::updateDt(){//calculate _dt
  //_dt=2.0/(_lambda*_rho*(double(_n)-1.0)*_Yhat);//
  _dt=2.0*double(_n)/(_lambda*_rho*double(_n-_Nc)*double(_n+_Nc-1)*_Yhat);
}


void DSMC::getnc(){//calculate nc of physical condensate
  _nc=double(_Nc)*_rho/((double)_n);
}


void DSMC::getNc(){//calculate Nc of simulated particles
  _Nc=0;
  for(int i=0;i<_n;i++){
    if(_pVec[i][3]<_pmin) _Nc++;
  }
}
