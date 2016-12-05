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

#define CONSIZE 2

// multiplication of two power series expansions = convolution of coeffcs
// this version is lobotomized to do only expansions up to order one.
// The t's are the expansion coefficients of a matrix
// the v's are the expansion coefficients of a vector.
// The result is stored back in the v's.

int convolve( complex<double> *t0, complex<double> *t1,
	       complex<double> *v0, complex<double> *v1)
{
  complex<double> y0[CONSIZE], y1[CONSIZE];

  y0[0] = t0[0]*v0[0] + t0[1]*v0[1];
  y0[1] = t0[2]*v0[0] + t0[3]*v0[1];

  y1[0] = t0[0]*v1[0] + t0[1]*v1[1];
  y1[1] = t0[2]*v1[0] + t0[3]*v1[1];

  y1[0] += t1[0]*v0[0] + t1[1]*v0[1];
  y1[1] += t1[2]*v0[0] + t1[3]*v0[1];

  v0[0] = y0[0];
  v0[1] = y0[1];
  v1[0] = y1[0];
  v1[1] = y1[1];

  return 0;
}

