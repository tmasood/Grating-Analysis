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

// initialize sheets. 

int initsheets(int **sheetsleft, int **sheetsright, complex<double> gam,
	       grating *gratptr, structure *epiptr)
{
  int i, k0, k;
  int order, ord2;
  int *shtsleft, *shtsright;

  double twopi, period;
  double widthl;
  double omega, wvl;
  double indx_r, indx_i;
  double loss;

  complex<double> kapleft, kapright;

  layer *layleftptr, *layrightptr;

  twopi = 2.0 * M_PI;
  // get grating period
  period = gratptr->getperiod();

  widthl = period/twopi;

  order = gratptr->getssspchrmncs();
  ord2 = (2 * order) + 1;

  shtsleft  = new int[ord2];
  shtsright = new int[ord2];

  layleftptr = epiptr->layerptr;
  layrightptr = epiptr->layerptr;
  wvl = epiptr->getwavelength();
  omega = (2 * M_PI)/wvl;
  indx_r = layleftptr->getlayerindex();
  loss = layleftptr->getlayerloss();
  // calculate imaginary part of index from loss
  indx_i = ((loss * wvl)/(4 * M_PI * 10000));
  kapleft = complex<double>(indx_r, indx_i);
  kapleft *= omega;

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

  for (i=0; i<ord2; i++)
    {
      shtsleft[i] = 1;
      shtsright[i] = 1;
    }

  if (real(gam) > 0.0)
    {
      k0 = (int)floor(-(real(kapleft) + imag(gam))*widthl);
      for (k=-k0; k>=-order; k--)
	{
	  shtsleft[k+order] = -1;
	}
      k0 = (int)floor(-(real(kapright) + imag(gam))*widthl);
      for (k=-k0; k>=-order; k--)
	{
	  shtsright[k+order] = -1;
	}
    }
  else
    {
      k0 = (int)floor((real(kapleft) - imag(gam))*widthl) + 1;
      for (k=k0; k<=order; k++)
	{
	  shtsleft[k+order] = -1;
	}
      k0 = (int)floor((real(kapright) - imag(gam))*widthl) + 1;
      for (k=k0; k<=order; k++)
	{
	  shtsright[k+order] = -1;
	}
    }

  *sheetsleft = shtsleft;
  *sheetsright = shtsright;
  return 0;
}
