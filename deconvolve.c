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
#include <cstdlib>

using namespace std;

// division of two power series expansions = deconvolution of coeffcs
// this version is lobotomized to do only expansions up to order one.
// the 0th components of v0 and v1 are in the denominator
// the 1st components of v0 and v1 are in the numerator
// res contains the expansion coefficients of the result

int deconvolve(complex<double> *v0, complex<double> *v1,
		complex<double> *res)
{
  complex<double> frac;

  if (abs(v0[0]) < 1e-10)
    {
      cout << "\n***deconvolve: division by zero\n" << endl;
      exit(1);
    }

  frac = 1.0/v0[0];
  res[0] = frac*v0[1];
  res[1] = (v1[1]*v0[0] - v1[0]*v0[1])*(frac*frac);

  return 0;
}

