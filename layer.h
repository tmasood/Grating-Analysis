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
#include <string.h>
#include <iostream>
#include <complex>

using namespace std;

#ifndef LAYER_INCLUDE
#define LAYER_INCLUDE

#define MAXCHAR 512

class layer
{
  double index;// Refractive index of layer
  double loss; // loss (1/cm)
  double lt;   // layer thickness
  double nlt;  // wavelength normalized layer thickness
  complex<double> kappa;
  complex<double> hz;
  complex<double> czh;
  complex<double> szh;
  char comments[MAXCHAR];

 public:
  layer *nextptr;
  layer *prevptr;

  layer();
  layer(double indx, double alpha, double th, char *layerinfo);
  ~layer();
  void setlayerkappa(complex<double> kap);
  void setlayerkth(double nt);
  void sethz(complex<double> ac);
  void setczh(complex<double> coszh);
  void setszh(complex<double> sinzh);
  double getlayerindex();
  double getlayerloss();
  double getlayerthickness();
  double getlayerkth();
  complex<double> gethz();
  complex<double> getczh();
  complex<double> getszh();
  complex<double> getlayerkappa();
};

#endif
