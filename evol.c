#include "dsmc.h"
#include <iostream>
using namespace std;

int main(){

  int  n=10000;
  //double sum=0,ep=1.0;
  double tMax=5.0e1,t=1.0;
  double tNext;//0.005039302255187421;		// time for save: corresponds to Francois' Delta t = 2000.   
  //n=3;
  double rho=0.2,pmin=0.1, m=0.5;
  DSMC expand(n,30,pmin,6.0,rho,t,m);

  tNext=expand.gettNext();
  do{
    if(expand._t>=tNext){
      expand.output();
      tNext=expand.gettNext();
      //cout << "tNext: " << tNext << endl;
     }
    expand.nextTime();
  }while(expand._t < tMax);
  expand.compTime();

}
