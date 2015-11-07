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
#include <fstream>
#include <sstream>
#include <complex>
#include "grating.h"
#include "gtoothpnl.h"
#include "layer.h"
#include "structure.h"

using namespace std;

// for a given point (x,z) in the interior domain, calculate the
// solution (=field value) using Green's representation theorem.
// The solution on the interfaces must be calculated by calcRhsInt()

complex<double> calcsolint( panel* pnls, domain *dmns, double x, double z)
{
  int k, j, i0, i1;
  domain *dmn;
  panel *pnl;
  double xC[2], alf, inner;
  Complex ptl[2], u, factor=0.25*imagI;

  xC[0] = x;  xC[1] = z;
  k = domainNr( x, z, &alf );
  if ( k==0 ) {
    printf("\n (%lf %lf) is not in the interior domain\n", x, z);
    exit(1);
  }
  factor *= 1.0/alf;
  
  for ( dmn=dmns, j=1; j<k; j++ ) dmn = dmn->next;

  for ( u=0.0, j=0; j<dmn->nPnls; j++ ) {
    pnl = &(pnls[dmn->pnls[j]]);
    i0 = pnl->idxU;
    i1 = pnl->idxV;
    calcp( pnl, xC, dmn->kappa, dmn->ornt[j], ptl);
    /* this is to cancel the self-term if there is one */
    inner = (x - pnl->xColl[0])*pnl->nrm[0] + (z - pnl->xColl[1])*pnl->nrm[1];
    if ( fabs(inner) < EPS ) ptl[1] = 0.0;

    if ( pnl->type <= 1 ) {
      u += dmn->ornt[j]*ptl[0]*Sol0[i1] - ptl[1]*Sol0[i0];
    }
    else if ( pnl->type == 2 ) {
      u -= lambda*( ptl[0]*Sol0[i1] + ptl[1]*Sol0[i0]);
    }
    else {
      u += ptl[0]*Sol0[i1] - ptl[1]*Rhs[i1] ;
    }
  }
  return factor*u;
} /* calcSolInt */
