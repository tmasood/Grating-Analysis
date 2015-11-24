//
//  Copyright 2015  Taha Masood, Johannes Tausch and Jerome Butler
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

#include <iostream>
#include <complex>

#include "structure.h"
#include "grating.h"

using namespace std;

// calculate the coefficients of the translation operators for given
// layer and a given mode.

int setuptransops(complex<double> *t0, complex<double> *t1,
		  int k, double widthlayer, complex<double> kap,
		  structure *epiptr, grating *gratptr)
{
  complex<double> ztmp, gkj, dgkj, ckj, skj;
  complex<double> gam;
  complex<double> twoPi;

  double period;

  twoPi = complex<double>(0.0,(2.0 * M_PI));

  // get grating period
  period = gratptr->getperiod();

  // get beta guess
  gam = epiptr->getbetaguess();

  ztmp = gam + twoPi*static_cast<double>(k)/period;
  gkj = sqrt(pow(ztmp,2.0) + pow(kap,2.0));
  dgkj = ztmp/gkj;
  ckj = cos(gkj*widthlayer);
  skj = sin(gkj*widthlayer);

  t0[0] = ckj;
  t0[1] = skj/gkj;
  t0[2] = -skj*gkj;
  t0[3] = t0[0];

  t1[0] = -skj*widthlayer*dgkj;
  t1[1] = (ckj*widthlayer - skj/gkj)*(dgkj/gkj);
  t1[2] = -(skj + ckj*widthlayer*gkj)*dgkj;
  t1[3] = t1[0];

  return 0;
}
