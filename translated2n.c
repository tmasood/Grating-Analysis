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
#include "grating.h"

using namespace std;

#define CONSIZE 2

int setuptransops(complex<double> *t0, complex<double> *t1,
		  int k, double widthlayer, complex<double> kap,
		  structure *epiptr, grating *gratptr);

int convolve( complex<double> *t0, complex<double> *t1,
	      complex<double> *v0, complex<double> *v1);

int deconvolve(complex<double> *v0, complex<double> *v1,
	       complex<double> *res);


// Translate the Dirichlet to Neumann operator from the infinite layer
// across layers. This can be done for the left as well as the right side.
// This process is described in the JCP paper, Sec 5, but is done here
// only for the expansion with linear terms.

// Parameters:
//    k           Fourier mode in z direction
//    g           on input D2N expansion of the inifinite layers
//                on return D2N expansion of the innermost layer
//    nJ          number of layers to translate
//    width       widths of each layer (positive on left, negative on right
//    kap         wavenumbers for each layer

int translated2n(int k, complex<double> *g, layer *layptr, int nj,
		 structure *epiptr, grating *gratptr)
{
  int j;
  complex<double> v0[CONSIZE], v1[CONSIZE];
  complex<double> t0[2*CONSIZE], t1[2*CONSIZE];
  complex<double> kap;
  double widthlayer;
  double omega, wvl;
  double indx_r, indx_i;
  double loss;

  v0[0] = 1;
  v0[1] = g[0];
  v1[0] = 0;
  v1[1] = g[1];

  layptr = layptr->nextptr;
  for (j = 1; j <= nj; j++)
    {
      wvl = epiptr->getwavelength();
      omega = (2 * M_PI)/wvl;
      indx_r = layptr->getlayerindex();
      loss = layptr->getlayerloss();
      // calculate imaginary part of index from loss
      indx_i = ((loss * wvl)/(4 * M_PI * 10000));
      kap = complex<double>(indx_r, indx_i);
      kap *= omega;
      widthlayer = layptr->getlayerthickness();
      setuptransops(t0, t1, k, widthlayer, kap, epiptr, gratptr);
      convolve(t0, t1, v0, v1);

      layptr = layptr->nextptr;
    }

  deconvolve(v0, v1, g);

  return 0;
}
