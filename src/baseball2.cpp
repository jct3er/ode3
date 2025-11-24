///
/// Starter template for second baseball problem
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
#include <TAxis.h>
#include <iostream>
#include <cstdio>
#include <cstdlib>

using namespace std;
const double g=32.152;
const double theta=0.01745;
const double B=0.00041;
const double omega=188.496;
const double m = 1.0;
double v0;
double phi;

double F(const vector<double> &y){
    return (0.0039 + 0.0058/(1+TMath::Exp((sqrt(y[1]*y[1] + y[3]*y[3] + y[5]*y[5]) -114.8294)/16.4042)))/2;
}

double f_ri(double x, const vector<double> &y, void *params=0){ 
  (void) x;   // prevent unused variable warning
  return y[1];
}

double f_vi(double x, const vector<double> &y, void *params){ 
  (void) x;   // prevent unused variable warning
  return -F(y)*sqrt(y[1]*y[1] + y[3]*y[3] + y[5]*y[5])*y[1] + B*omega*(y[5]*TMath::Sin(phi)-y[3]*TMath::Cos(phi));
}

double f_rj(double x, const vector<double> &y, void *params=0){ 
  (void) x;   // prevent unused variable warning
  return y[3];
}

double f_vj(double x, const vector<double> &y, void *params){ 
  (void) x;
  return -F(y)*sqrt(y[1]*y[1] + y[3]*y[3] + y[5]*y[5])*y[3] + B*omega*y[1]*TMath::Cos(phi);
}

double f_rk(double x, const vector<double> &y, void *params=0){ 
  (void) x;   // prevent unused variable warning
  return y[5];
}

double f_vk(double x, const vector<double> &y, void *params){ 
  (void) x;
  return -g-F(y)*sqrt(y[1]*y[1] + y[3]*y[3] + y[5]*y[5])*y[5] - B*omega*y[1]*TMath::Sin(phi);
}

double f_stop(double x, const vector<double> &y, void *params){
  (void) x;
  if (y[0]>60) return 1;  // stop calulation if the current step takes x past home plate
  return 0;  // continue calculation
}


int main(int argc, char **argv){

  // we have 6 initial conditions for this problem
  // y[0] = y[2] = y[4] = 0;  // init x,y,z
  // y[1] = v0*cos(theta0);   // vx  "x is line towards the plate
  // y[3] = 0;                // vy  "y" is measured as left/right divergence from line to plate
  // y[5] = v0*sin(theta0);   // vz  "z" is vertival measure
  vector<double> y0(6);

  bool showPlot=false;
  // pitches
  // slider ip=0
  // curve ip=1
  // screwball ip=2
  // fast ip=3
  int ip=1;    // default pitch
  int c;
  while ((c = getopt (argc, argv, "p:n")) != -1)
    switch (c) {
    case 'p':
      ip = atoi(optarg);
      break;
    case 'n':
      showPlot=true;
      break;
    }

  TString title;
  if (ip==0){
    cout << "Setting up initial conditions for slider" << endl;
    //SetupSlider(y0);
      phi = 0;
      v0 = 124.667;
  }
  else if (ip==1){
    cout << "Setting up initial conditions for curveball" << endl;
    //SetupCurve(y0);
      phi = 0.7854;
      v0 = 124.667;
  }
  else if (ip==2){
    cout << "Setting up initial conditions for screwball" << endl;
    //SetupScrewball(y0);
      phi = 2.3562;
      v0 = 124.667;
  }
  else {
    cout << "Setting up initial conditions for fastball" << endl;
    //SetupFastball(y0);
      phi = 3.927;
      v0 = 139.33;
  }

  TApplication theApp("App", &argc, argv); // init ROOT App for displays

  double xend=60;   // feet
  double yend=0;    // tbd
  double zend=0;    // tbd
  double vxend=0;
  double vyend=0;
  double vzend=0;

  // write code here

  y0[0]=0;
  y0[2]=0;
  y0[4]=0;
  y0[3]=0;
  y0[1]=v0*TMath::Cos(theta);
  y0[5]=v0*TMath::Sin(theta);

  vector<pfunc_t> v_fun(6);   // 4 element vector of function pointers
  v_fun[0]=f_ri;
  v_fun[1]=f_vi;
  v_fun[2]=f_rj;
  v_fun[3]=f_vj;
  v_fun[4]=f_rk;
  v_fun[5]=f_vk;

  double t=0;
  double h=0.0001;  // step size
  vector<TGraph> y = RK4SolveN(v_fun, y0, h, t, 0, f_stop, 10000);

  int nPoints = y[0].GetN();
  y[1].GetPoint(nPoints-1, t, vxend);
  y[2].GetPoint(nPoints-1, t, yend);
  y[3].GetPoint(nPoints-1, t, vyend);
  y[4].GetPoint(nPoints-1, t, zend);
  y[5].GetPoint(nPoints-1, t, vzend);

  TGraph xy;
  TGraph xz;

  double xpoint;
  double ypoint;
  double zpoint;

  for (int i=0;i<nPoints;i++){
    y[0].GetPoint(i, t, xpoint);
    y[2].GetPoint(i, t, ypoint);
    y[4].GetPoint(i, t, zpoint);
    xy.SetPoint(i,xpoint,ypoint);
    xz.SetPoint(i,xpoint,zpoint);
  }


  TCanvas* tc = new TCanvas();
  // xy.SetTitle("Fastball;X (ft);Z (ft) / Y (ft)");
  // xy.GetYaxis()->SetRangeUser(-4, 2);
  // xy.SetLineStyle(kDotted);
  // xy.SetLineWidth(1);
  // xy.Draw();
  // xz.Draw("SAME");
  // tc->Draw();
  // tc->Update();
  // tc->Print("../pitches_fastball.pdf", "Title: Fastball");

  if (ip==0){
      xy.SetTitle("Slider;X (ft);Z (ft) / Y (ft)");
      xy.GetYaxis()->SetRangeUser(-4, 2);
      xy.SetLineStyle(kDotted);
      xy.SetLineWidth(1);
      xy.Draw();
      xz.Draw("SAME");
      tc->Draw();
      tc->Update();
      tc->Print("../pitches_slider.pdf", "Title: Slider");
  }
  else if (ip==1){
      xy.SetTitle("Curveball;X (ft);Z (ft) / Y (ft)");
      xy.GetYaxis()->SetRangeUser(-4, 2);
      xy.SetLineStyle(kDotted);
      xy.SetLineWidth(1);
      xy.Draw();
      xz.Draw("SAME");
      tc->Draw();
      tc->Update();
      tc->Print("../pitches_curveball.pdf", "Title: Curveball");
  }
  else if (ip==2){
      xy.SetTitle("Screwball;X (ft);Z (ft) / Y (ft)");
      xy.GetYaxis()->SetRangeUser(-4, 2);
      xy.SetLineStyle(kDotted);
      xy.SetLineWidth(1);
      xy.Draw();
      xz.Draw("SAME");
      tc->Draw();
      tc->Update();
      tc->Print("../pitches_screwball.pdf", "Title: Screwball");
  }
  else {
      xy.SetTitle("Fastball;X (ft);Z (ft) / Y (ft)");
      xy.GetYaxis()->SetRangeUser(-4, 2);
      xy.SetLineStyle(kDotted);
      xy.SetLineWidth(1);
      xy.Draw();
      xz.Draw("SAME");
      tc->Draw();
      tc->Update();
      tc->Print("../pitches_fastball.pdf", "Title: Fastball");
  }


  printf("Time to pass plate:  (%lf)\n",h*nPoints);

  // to compare to the plots in Fitzpatrick, output your results in **feet**
  // do not change these lines
  printf("********************************\n");
  printf("Coordinates when x=60 feet\n");
  printf("(x,y,z) = (%lf,%lf,%lf)\n",xend,yend,zend);
  printf("(vx,vy,vz) = (%lf,%lf,%lf)\n",vxend,vyend,vzend);
  printf("********************************\n");

  // plot the trajectory.  See Fitzpatrick for plot details
  if (showPlot){      
    cout << "Press ^c to exit" << endl;
    theApp.SetIdleTimer(30,".q");  // set up a failsafe timer to end the program  
    theApp.Run();
  }
  
  return 0;
}

