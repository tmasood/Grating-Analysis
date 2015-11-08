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
#include "gtoothdmn.h"
#include "grating.h"

using namespace std;

// Calculate the collocation points for the boundary of each domain.
// Note that collocpoints are also part of the panels structure.
// Returns the total number of points = number of rows/cols in the matrix

int getcollocpts(grating *gratptr)
{
  gtdomain *dmn;
  int npts, pnls, nptsall;

  nptsall = 0;
  npts = 0;
  for (dmn = gratptr->gtdmnptr; dmn != NULL; dmn = dmn->nextptr)
    {
      pnls = dmn->getnpanels();
      npts = pnls;
      dmn->setnpts(npts);
      nptsall += npts;
    }

  return nptsall;
}
