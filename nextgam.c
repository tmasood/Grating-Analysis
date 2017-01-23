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

#include "structure.h"
#include "layer.h"
#include "grating.h"

using namespace std;

// return the infinity norm of the matrix stored columnwise

complex<double> nextgam(int order, complex<double> *z, complex<double> gam,
			int *sheetsleft, int *sheetsright, 
			structure *epiptr, grating *gratptr)
{
  int i, k, n;

  double normal, t0, t1, kap;
  double omega, wvl;
  double indx_r, indx_i;
  double loss;
  double period;

  complex<double> twopi;
  complex<double> gamnew, newtonstep, gamk;
  complex<double> kapleft, kapright;

  layer *layleftptr, *layrightptr;

  n = (4*order)+2;
  normal = 0;
  twopi = complex<double>(0.0,(2.0 * M_PI));

  period = gratptr->getperiod();

  // Newton update is the reciprocal of the largest eigenvalue
  for (i=0; i<n; i++)
    {
      if (abs(z[i]) > normal)
	{
	  newtonstep = 1.0 / z[i];
	  gamnew = gam + newtonstep;
	  normal = abs(z[i]);
	}
    }

  if (fabs(real(newtonstep)) > 1e-15)
    {
      // change sheets, if necessary, for the left
      layleftptr = epiptr->layerptr;
      wvl = epiptr->getwavelength();
      omega = (2 * M_PI)/wvl;
      indx_r = layleftptr->getlayerindex();
      loss = layleftptr->getlayerloss();
      // calculate imaginary part of index from loss
      indx_i = ((loss * wvl)/(4 * M_PI * 10000));
      kapleft = complex<double>(indx_r, indx_i);
      kapleft *= omega;
      kap = abs(kapleft);

    for (i=0, k=-order; k<=order; k++, i++)
      {
	gamk = gam + ((twopi*static_cast<double>(k))/period);
	t0 = -(real(gamk))/(real(newtonstep));
	if (t0 >= 0.0 && t0 <= 1.0)
	  {
	    t1 = imag(gamk) + t0*imag(newtonstep);
	    if (fabs(t1) >  kap)
	      {
		sheetsleft[i] *= -1;
	      }
	  }
      }
    
    // change sheets, if necessary, for the right
    layrightptr = epiptr->layerptr;

    while (layrightptr->nextptr != NULL)
      {
	layrightptr = layrightptr->nextptr;
      }
    indx_r = layrightptr->getlayerindex();
    loss = layrightptr->getlayerloss();
    // calculate imaginary part of index from loss
    indx_i = ((loss * wvl)/(4 * M_PI * 10000));
    kapright = complex<double>(indx_r, indx_i);
    kapright *= omega;
    kap = abs(kapright);
    for (i=0, k=-order; k<=order; k++, i++)
      {
	gamk = gam + (twopi*static_cast<double>(k))/period;
	t0 = - (real(gamk))/(real(newtonstep));
	if (t0 >= 0.0 && t0 <= 1.0)
	  {
	    t1 = imag(gamk) + t0*imag(newtonstep);

	    if ( fabs(t1) >  kap)
	      {
		sheetsright[i] *= -1;
	      }
	  }
      }
    }
  return gamnew;
}
