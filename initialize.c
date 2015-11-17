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

using namespace std;

int initialize(structure *epiptr)
{
  layer *layptr, *headlayptr;
  complex<double> c_index, c2_index;
  complex<double> zguess, z2guess, betaguess;
  complex<double> arg, ac;
  complex<double> sinzh, coszh;

  double loss;
  double index_r;
  double index_i;
  double wavelength;
  double guess_r;
  double guess_i;
  double th;
  double ko;
  double nlt;
  int numlayers;

  wavelength = epiptr->getwavelength();
  ko = (2 * M_PI)/wavelength;
  epiptr->setko(ko);
  // guess effective index squared
  guess_r = epiptr->getguessr();
  guess_i = epiptr->getguessi();
  zguess = complex<double>(guess_r, guess_i);
  z2guess = pow(zguess,2.0);
  epiptr->setz2guess(z2guess);
  betaguess = ko * zguess;
  epiptr->setbetaguess(betaguess);

  layptr = epiptr->layerptr;
  headlayptr = epiptr->layerptr;
  numlayers = 0;

  while (layptr != NULL)
    {
      // kappa (layer index squared)
      index_r = layptr->getlayerindex();
      loss = layptr->getlayerloss();
      // calculate imaginary part of index from loss
      index_i = ((loss * wavelength)/(4 * M_PI * 10000));
      c_index = complex<double>(index_r, index_i);
      c2_index = pow(c_index,2.0);
      // cout << c2_index << endl;
      // store kappa
      layptr->setlayerkappa(c2_index);
      // normalized layer thickness
      th = layptr->getlayerthickness();
      nlt = ko * th;
      layptr->setlayerkth(nlt);
      // hz
      if ((headlayptr == layptr) || (layptr->nextptr == NULL))
	{
	  // the first layer or the last layer hz
	  ac = sqrt(-1.0 * (z2guess + c2_index));
	  layptr->sethz(ac);
	}
      else
	{
	  ac = sqrt(1.0 * (z2guess + c2_index));
	  layptr->sethz(ac);
	  arg = nlt * ac;
	  coszh = cos(arg);
	  layptr->setczh(coszh);
	  sinzh = sin(arg);
	  layptr->setszh(sinzh);
	}

      numlayers++;
      layptr = layptr->nextptr;
    }

  epiptr->settotallayers(numlayers);

  return 0;
}
