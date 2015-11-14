//
//  Copyright 2015 Quantum Designs LLC, Taha Masood, Johannes Tausch
//  and Jerome Butler
//
//  Permission to use, copy, and distribute this software and its
//  documentation for any purpose with or without fee is hereby granted,
//  provided that the above copyright notice appear in all copies and
//  that both the copyright notice and this permission notice appear
//  in supporting documentation.
//
//  This software is provided "as is" without express or implied warranty
//  to the extent permitted by applicable law.
//
#include <complex>
#include "gtoothpnl.h"
#include "gtoothdmn.h"

#ifndef GRATING_INCLUDE
#define GRATING_INCLUDE

class grating
{
  int grtlayer;
  int fillayer;
  double period;
  double dutycycle;
  double dx;  // sub panel length
  int ssspchrmncs;  // single sided space harmonics
  int cvertpnts; // common vertical points
  int totalpnls; // total refined panels
  string gshape; // grating shape

 public:
  gtpanel *gtpnlptr; // pointer to the grating tooth panels
  gtpanel *gtrefpnlptr; // pointer to the refined grating tooth panels
  gtdomain *gtdmnptr; // pointer to the grating tooth domains
  complex<double> *gzk;
  
  grating();
  ~grating();
  void setgllayer(int gl);
  void setfllayer(int fl);
  void setperiod(double prd);
  void setdc(double dc);
  void setdx(double pl);
  void setssspchrmncs(int sssh);
  void setcvertpnts(int cvpts);
  void setgratingshape(string gs);
  void settotalpnls(int tpnls);

  int getgl(); // get the grating layer number
  int getfl(); // get the fill layer number
  string getgratingshape(); // get the grating shape
  double getperiod();
  double getdc();
  double getdx();
  int getssspchrmncs();
  int gettotalpnls();
};

#endif
