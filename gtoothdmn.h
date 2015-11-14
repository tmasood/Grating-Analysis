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
// gtdomain - Grating tooth domains.
//
#include <iostream>
#include <complex>

using namespace std;

#ifndef GTDOMAIN_INCLUDE
#define GTDOMAIN_INCLUDE

#define sgn(a) (((a) >= 0 ) ? 1 : (-1))

class gtdomain
{
  complex<double> kap;// layer domain index
  int npnls; // number of panels that border the domain
  int npts;

 public:
  int *indp; // indices of panels on the boundary
  int *ornt; // orientation of panel normals
  gtdomain *nextptr;

  gtdomain();
  ~gtdomain();
  void setkap(complex<double> kappa);
  void setnpanels(int npanels);
  void setnpts(int pts);

  complex<double> getkap();
  int getnpanels();
  int getnpts();
};

#endif
