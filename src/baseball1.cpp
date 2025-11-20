///
/// Starter template for first baseball problem
/// Solve for the initial speed of the pitch given the initial parameters
/// xend : distance to home plate [18.5] m
/// z0 : height of release of ball [1.4] m
/// theta0 : angle of release above horizontal [1] degree
///
///  Do not change the interface for running the program
///  Fill in the value of vPitch in the print statement with your solution
///  at the end of main()
///

#include "RKn.hpp"
#include "TROOT.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGClient.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>

using namespace std;

double xend;

struct Params {
  double g;   // acceleration [m/s^2]
  double m;   // mass of object [kg], nb proj. In vacuum funcs do not depend on the mass
  double d;   // m diameter of ball
  double b;   // b,c params for air resistance
  double c;
};



double f_ri(double x, const vector<double> &y, void *params=0){ 
  (void) x;   // prevent unused variable warning
  return y[1];
}

double f_vi(double x, const vector<double> &y, void *params){ 
  (void) x;   // prevent unused variable warning
  Params* p = (Params*) params;
  return -(p->b*p->d*y[1] + p->c*p->d*p->d*sqrt(y[1]*y[1] + y[3]*y[3])*y[1])/p->m;
}

double f_rj(double x, const vector<double> &y, void *params=0){ 
  (void) x;   // prevent unused variable warning
  return y[3];
}

double f_vj(double x, const vector<double> &y, void *params){ 
  (void) x;   // prevent unused variable warning
  Params* p = (Params*) params;
  return -(p->b*p->d*y[3] + p->c*p->d*p->d*sqrt(y[1]*y[1] + y[3]*y[3])*y[3])/p->m - p->g;
}

double f_stop(double x, const vector<double> &y, void *params){
  (void) x;
  if (y[0]>xend) return 1;  // stop calulation if the current step takes x past home plate
  return 0;  // continue calculation
}


int main(int argc, char **argv){

  // examples of parameters
  Params pars;
  pars.g=9.81;
  pars.m=0.145;    
  pars.d=0.075;   
  pars.b=1.6e-4;  
  pars.c=0.25;
  void *p_par = (void*) &pars;

  xend=18.5;       // meters to plate
  double z0=1.4;             // height of release [m]
  double theta0=1;         // angle of velocity at release (degrees)
                                      // convert to radians before using!
  bool showPlot=false;    // keep this flag false by default
  
  // allow changing the parameters from the command line
  int c;
  while ((c = getopt (argc, argv, "x:z:t:p")) != -1)
    switch (c) {
    case 'x':
      xend = atof(optarg);
      break;
    case 'z':
      z0 = atof(optarg);
      break;
    case 't':
      theta0 = atof(optarg);
      break;
    case 'p':
      showPlot=true;
      break;
    case '?':
      fprintf (stderr, "Unknown option `%c'.\n", optopt);
    }
  TApplication theApp("App", &argc, argv); // init ROOT App for displays


  double vPitch = 40;   // m/s of pitch needed to land in strike zone at 0.9 meters
  // write code to solve for vPitch here

  vector<double> y0(4);
  // y0[0] = 0;
  // y0[2] = z0;
  // y0[1] = vPitch*TMath::Cos(theta0*3.1415/180);
  // y0[1] = vPitch*TMath::Sin(theta0*3.1415/180);

  vector<pfunc_t> v_fun(4);   // 4 element vector of function pointers
  v_fun[0]=f_ri;
  v_fun[1]=f_vi;
  v_fun[2]=f_rj;
  v_fun[3]=f_vj;

  double t;
  double h=0.0001;  // step size
  vector<TGraph> y;

  double zend=0;
  int nPoints;
  // y[2].GetPoint(nPoints-1, t, zend);


  while (zend<0.9){
    vPitch += .01;
    y0[0] = 0;
    y0[2] = z0;
    y0[1] = vPitch*TMath::Cos(theta0*3.1415/180);
    y0[3] = vPitch*TMath::Sin(theta0*3.1415/180);

    t=0;
    y = RK4SolveN(v_fun, y0, h, t, p_par, f_stop, 10000);
    nPoints = y[0].GetN();
    y[2].GetPoint(nPoints-1, t, zend);
    //printf("vPitch, zend = (%lf,%lf)\n", vPitch, zend);
  }




  // do not change these lines
  printf("********************************\n");
  printf("(xend,z0,theta0) = (%lf,%lf,%lf)\n",xend,z0,theta0);
  printf("v_pitch = %lf m/s\n",vPitch);
  printf("********************************\n");

  if (showPlot){
    cout << "Press ^c to exit" << endl;
    theApp.SetIdleTimer(30,".q");  // set up a failsafe timer to end the program  
    theApp.Run();
  }
  
  return 0;
}

