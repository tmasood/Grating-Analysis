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

// define domain for sawtooth grating
int defdmnst(structure *epiptr, grating *gratptr)
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

  // top domain consists of fill region
  // point to the grating fill layer
  for (i = 0; i < fl; i++)
    {
       glayptr = glayptr->nextptr;
    }

 // get the grating fill layer index
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
  dmnptr->setnpanels(3);
  dmnptr->indp = new int[3];
  // indices of panels on boundary
  dmnptr->indp[0] = 7;
  // panels whose normal points into the domain
  dmnptr->indp[1] = 1;
  // other 2 panels
  dmnptr->indp[2] = 3;

  // middle interior domain consists of grating region
  dmnptr->nextptr = new gtdomain;
  dmnptr = dmnptr->nextptr;

  glayptr = epiptr->layerptr;

  // point to the grating layer
  for (i = 0; i < gl; i++)
    {
       glayptr = glayptr->nextptr;
    }

 // get the grating layer index
  indx_r = glayptr->getlayerindex();
  loss = glayptr->getlayerloss();
  // calculate imaginary part of index from loss
  indx_i = ((loss * wvl)/(4 * M_PI * 10000));
  kappa = complex<double>(indx_r, indx_i);
  kappa *= omega;
  dmnptr->setkap(kappa);

  // number of panels that border the domain
  dmnptr->setnpanels(4);
  dmnptr->indp = new int[4];
  // indices of panels on boundary
  dmnptr->indp[0] = -2;
  // panels whose normal point into the domain
  dmnptr->indp[1] = -3;
  // panel 2 is common between both domains. Minus sign because
  // normal points in different direction
  dmnptr->indp[2] = 5;
  dmnptr->indp[3] = 4;

  // bottom domain consists of fill region

  // point to the grating fill layer
  for (i = 0; i < fl; i++)
    {
       glayptr = glayptr->nextptr;
    }

 // get the grating fill layer index
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
  dmnptr->setnpanels(3);
  dmnptr->indp = new int[3];
  // indices of panels on boundary
  dmnptr->indp[0] = 0;
  // panels whose normal point into the domain
  dmnptr->indp[1] = 6;
  // panel 2 is common between both domains. Minus sign because
  // normal points in different direction
  dmnptr->indp[2] = 2;

  // end of domain list
  dmnptr->nextptr = NULL;

  return 0;
}
