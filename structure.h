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
#include <iomanip>
#include "layer.h"
#include "gtoothpnl.h"
#include "grating.h"

using namespace std;

#ifndef STRUCTURE_INCLUDE
#define STRUCTURE_INCLUDE

class structure
{
  double wavelength;
  double zguessi;
  double zguessr;
  double ko;
  int totallayers;
  complex<double> z2guess;
  complex<double> betaguess;
  string structtype;

 public:
  layer  *layerptr; // pointer to the current layer structure
  grating *gratptr; // pointer to the grating structure
  structure *nextptr; // pointer to next structure

  structure();
  ~structure();
  void setwavelength(double wvl);
  void setguess(double zgi, double zgr);
  void setz2guess(complex<double> z2g);
  void setbetaguess(complex<double> betag);
  void setko(double wv);
  void settotallayers(int numlayer);
  void setstructtype(string st);
  double getwavelength();
  double getko();
  double getguessr();
  double getguessi();
  int gettotallayers();
  string getstructtype();
  complex<double> getbetaguess();
};

#endif
