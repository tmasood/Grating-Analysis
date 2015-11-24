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
#include "structure.h"

// default constructor
structure::structure()
{
  wavelength = 1.55;
  zguessi = 0.0;
  zguessr = 0.0;
  z2guess = complex<double>(0.0, 0.0);
  betaguess = complex<double>(0.0, 0.0);
  totallayers = 0;
  structtype = "WG";
  layerptr = NULL;
  gratptr = NULL;
  nextptr = NULL;
}

// default destructor
structure::~structure()
{
  delete [] layerptr;
  delete [] gratptr;
  delete [] nextptr;
}

void structure::setwavelength(double wvl)
{
  wavelength = wvl;
}

void structure::setguess(double zgi, double zgr)
{
  zguessi = zgi;
  zguessr = zgr;
}

void structure::setz2guess(complex<double> z2g)
{
  z2guess = z2g;
}

void structure::setbetaguess(complex<double> betag)
{
  betaguess = betag;
}

void structure::setko(double wv)
{
  ko = wv;
}

void structure::settotallayers(int numlayers)
{
  totallayers = numlayers;
}

void structure::setstructtype(string st)
{
  structtype = st;
}

double structure::getwavelength()
{
  return wavelength;
}

double structure::getko()
{
  return ko;
}


double structure::getguessr()
{
  return zguessr;
}

double structure::getguessi()
{
  return zguessi;
}

int structure::gettotallayers()
{
  return totallayers;
}

string structure::getstructtype()
{
  return structtype;
}

complex<double> structure::getbetaguess()
{
  return (betaguess);
}
