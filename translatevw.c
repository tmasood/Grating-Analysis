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

// layers are coupled via the relation described in translatevw

int  translatevw(int order, complex<double> *v, complex<double> *w,
		 complex<double> *gkb, double width)
{
  int k, ord2;
  complex<double> v1, w1;

  ord2=2*(order+1);
  for (k=0; k < ord2; k++)
    {
      v1 = cos(gkb[k]*width) * v[k];
      v1 += (sin(gkb[k]*width) / gkb[k]) * w[k];
      w1 = cos(gkb[k]*width) * w[k];
      w1 -= gkb[k] * sin(gkb[k]*width) * v[k];

      // state variables for the next adjacent layer
      v[k] = v1;
      w[k] = w1;
    }

  return 0;
} 
