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
#include "grating.h"

// default constructor
grating::grating()
{
  grtlayer = 0;
  fillayer = 0;
  period = 0.0;
  dutycycle = 0.0;
  dx = 0.0;
  ssspchrmncs = 0;
  cvertpnts = 0;
  gtpnlptr = NULL;
  gtrefpnlptr = NULL;
  gtdmnptr = NULL;
  gzk = NULL;
}

// default destructor
grating::~grating()
{
  delete [] gtpnlptr;
  delete [] gtdmnptr;
}

void grating::setgllayer(int gl)
{
  grtlayer = gl;
}

void grating::setfllayer(int fl)
{
  fillayer = fl;
}

void grating::setperiod(double prd)
{
  period = prd;
}

void grating::setdc(double dc)
{
  dutycycle = dc;
}

void grating::setdx(double pl)
{
  dx = pl;
}

void grating::setssspchrmncs(int sssh)
{
  ssspchrmncs = sssh;
}

void grating::setcvertpnts(int cvpts)
{
  cvertpnts = cvpts;
}

void grating::settotalpnls(int tpnls)
{
  totalpnls = tpnls;
}

void grating::setgratingshape(string gs)
{
  gshape = gs;
}

int grating::getgl()
{
  return grtlayer;
}

int grating::getfl()
{
  return fillayer;
}

string grating::getgratingshape()
{
  return gshape;
}

double grating::getperiod()
{
  return period;
}

double grating::getdc()
{
  return dutycycle;
}

double grating::getdx()
{
  return dx;
}

int grating::getssspchrmncs()
{
  return ssspchrmncs;
}

int grating::gettotalpnls()
{
  return totalpnls;
}
