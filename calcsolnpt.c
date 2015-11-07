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
#include <iostream>
#include <complex>

using namespace std;

// evaluate the solution which is given by
//   u(x,0) = \sum_k v_k cos( g_k h) + w_k sin(g_k h)
// this can be used for positive and negative h

complex<double> calcsolnpt(int order, complex<double> *v, complex<double> *w,
			   complex<double> *gkb, double h)
{
  complex<double> fval;
  int ord2;
  int k;

  fval = 0.0;
  ord2 = 2*order+1;

  // from equation 24 of Dr Johannes Tausch Paper
  // u(x,0) = Sum_over_n[phi_n * cos(kxn(x - xj)) +
  // chi_n * (u/kxn) * sin(kxn(x - xj))]

  // thus for each harmonic kxn = gkb
  // (x - xj) = h, v = phi, and w = chi
  for (k = 0; k < ord2; k++)
    {
      fval += v[k]*cos(gkb[k]*h);
      fval += w[k]*sin(gkb[k]*h)/gkb[k];
    }
  return fval;
}
