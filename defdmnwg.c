//
//  Copyright 2016 Quantum Designs LLC, Taha Masood, Johannes Tausch
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
#include <cstdlib>

#include "gtoothdmn.h"
#include "grating.h"
#include "structure.h"

using namespace std;

int defdmnwg(structure *epiptr, grating *gratptr)
{
  layer *glayptr;
  gtdomain *dmnptr;

  int gl, fl;
  int i;

  complex<double> kappa;

  double indx_r, indx_i;
  double loss;
  double omega, wvl;

  dmnptr = new gtdomain;
  gratptr->gtdmnptr = dmnptr;

  glayptr = epiptr->layerptr;
  // get grating layer
  gl = gratptr->getgl();
  fl = gratptr->getfl();

  if ((fl > (gl + 1)) || (fl < (gl - 1)))
    {
      cout << endl;
      cout << "fill layer must be adjacent to grating layer " << endl;
      cout << "...... exiting ......" << endl;
      exit(1);
    }

  // point to the grating layer
  for (i = 0; i < gl; i++)
    {
       glayptr = glayptr->nextptr;
    }

 // get the grating layer index
  indx_r = glayptr->getlayerindex();
  loss = glayptr->getlayerloss();
  wvl = epiptr->getwavelength();
  // calculate imaginary part of index from loss
  indx_i = ((loss * wvl)/(4 * M_PI * 10000));
  kappa = complex<double>(indx_r, indx_i);
  omega = (2 * M_PI)/wvl;
  kappa *= omega;
  dmnptr->setkap(kappa);

  // number of panels that border the domain
  dmnptr->setnpanels(4);
  dmnptr->indp = new int[4];
  // indices of panels on boundary
  dmnptr->indp[0] = 0;
  // panels whose normal point into the domain
  dmnptr->indp[1] = 1;
  // panel 2 is common between both domains. Minus sign because
  // normal points in different direction
  dmnptr->indp[2] = 2;
  dmnptr->indp[3] = 3;

  // end of domain list
  dmnptr->nextptr = NULL;

  return 0;
}
