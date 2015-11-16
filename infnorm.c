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

// return the infinity norm of the matrix stored columnwise

double infnorm(int n, complex<double> *a)
{
  int i, j;
  double maxval=0, tmp;

  for (i=0; i<n; i++)
    {
      for (tmp=0,j=0; j<n; j++)
	{
	  tmp += abs(a[i + n*j]);
	}
      if (tmp > maxval)
	{
	  maxval = tmp;
	}
    }
  return maxval;
}
