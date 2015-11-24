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

#define CONSTSIZE 2

int translated2n(int k, complex<double> *g, layer *layptr, int gl,
                 structure *epiptr, grating *gratptr);
int translated2nr(int k, complex<double> *g, layer *layptr, int gl,
                 structure *epiptr, grating *gratptr);

// subtract exterior D2N map from the moment matrices

int subext(complex<double> **mt, complex<double> **d2nleft,
	    complex<double> **d2nright, structure *epiptr,
	    grating *gratptr, int *sheetsleft, int *sheetsright)
{
  int order, ord2, ord4;
  int k, k1, k2, m;
  complex<double> gkbl, gkbr, I;
  complex<double> kapleft2, kapright2, tmp;
  complex<double> kapleft, kapright;
  complex<double> gleft[CONSTSIZE], gright[CONSTSIZE];
  complex<double> gam;
  complex<double> *extd2nlft, *extd2nrght;
  complex<double> twopi;

  double omega, wvl;
  double indx_r, indx_i;
  double loss;
  int gl, rgl;
  double period;

  layer *layleftptr, *layrightptr;

  extd2nlft = NULL;
  extd2nrght = NULL;

  order = gratptr->getssspchrmncs();
  ord2 = (2 * order) + 1;
  ord4 = (2 * ord2);
  real(I) = 0.0;
  imag(I) = 1.0;
  layleftptr = epiptr->layerptr;
  layrightptr = epiptr->layerptr;
  while (layrightptr->nextptr != NULL)
    {
      layrightptr = layrightptr->nextptr;
    }

  wvl = epiptr->getwavelength();
  omega = (2 * M_PI)/wvl;
  indx_r = layleftptr->getlayerindex();
  loss = layleftptr->getlayerloss();
  // calculate imaginary part of index from loss
  indx_i = ((loss * wvl)/(4 * M_PI * 10000));
  kapleft = complex<double>(indx_r, indx_i);
  kapleft *= omega;
  kapleft2 = (kapleft * kapleft);

  // get grating layer
  gl = gratptr->getgl();
  rgl = (epiptr->gettotallayers() - gl - 2);

  indx_r = layrightptr->getlayerindex();
  loss = layrightptr->getlayerloss();
  // calculate imaginary part of index from loss
  indx_i = ((loss * wvl)/(4 * M_PI * 10000));
  kapright = complex<double>(indx_r, indx_i);
  kapright *= omega;
  kapright2 = (kapright * kapright);

  if (extd2nlft == NULL)
    {
      extd2nlft = new complex<double>[ord2];
      extd2nrght = new complex<double>[ord2];
    }

  // set up the coefficients
  for (k = -order, k1 = 0, k2 = ord2; k <= order; k++, k1++, k2++)
    {
      // D2N and its derivative for infinite layers
      // tmp = kzn = (2*pi*n / period) + gamma
      // n is the order of harmonic.
      // gkbl and gkbr are kxn = sqrt(kzn**2 + u_j*epsilon_j*(omega**2))
      period = gratptr->getperiod();
      gam = epiptr->getbetaguess();
      twopi = complex<double>(0.0,(2.0 * M_PI));
      tmp = gam + (twopi*static_cast<double>(k)/period);
      gkbl = double(sheetsleft[k1])*sqrt(pow(tmp,2.0) + kapleft2);
      gkbr = double(sheetsright[k1])*sqrt(pow(tmp,2.0) + kapright2);

      gleft[0] = -I*gkbl;
      gleft[1] = -I*tmp/gkbl; // -1 because on the left x -> (-x)
      gright[0] = I*gkbr;
      gright[1] = I*tmp/gkbr;

      if (gl > 0)
        {
          translated2n(k, gleft, layleftptr, gl, epiptr, gratptr);
        }

      if (rgl > 0)
        {
          translated2nr(k, gright, layrightptr, rgl, epiptr, gratptr);
        }

      // subtract expansions
      for (m = 0; m <= 1; m++)
	{
	  mt[m][k1*ord4+k1] += gleft[m];  // +, because of normal derivative
	  mt[m][k2*ord4+k2] -= gright[m];
	}

      extd2nlft[k1] = -gleft[0];   // save D2N to calculate the solution
      extd2nrght[k1] = gright[0];
    }
  *d2nright = extd2nrght;
  *d2nleft = extd2nlft;

  return 0;
}
