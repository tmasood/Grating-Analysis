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
#include <cstdlib>

#include "grating.h"
#include "gtoothpnl.h"
#include "layer.h"
#include "structure.h"

using namespace std;

int domainnr(double x, double z, double *alf, structure *epiptr,
	     grating *gratptr);
int calcp(gtpanel *p, double *x, complex<double> kappa, int ornt,
           complex<double> *ptl);

extern double EPS;

// for a given point (x,z) in the interior domain, calculate the
// solution (=field value) using Green's representation theorem.
// The solution on the interfaces must be calculated by calcRhsInt()

complex<double> calcsolint( gtpanel *pnls, gtdomain *dmns, double x,
			    double z, complex<double> lambda,
			    complex<double> *rhs, complex<double> **sol0t,
			    structure *epiptr, grating *gratptr)
{
  int k, j, i0, i1;
  int pnlidxu, pnlidxv;
  int npnls;
  int type;

  gtdomain *dmn;
  gtpanel *pnl;
  double xc[2], alf, inner;
  double *xcoll, *nrm;

  complex<double> I;
  complex<double> ptl[2], u, factor;
  complex<double> kappa;
  complex<double> *sol0;

  I = complex<double>(0.0, 1.0);

  sol0 = *sol0t;

  factor=0.25*I;
  xc[0] = x;  xc[1] = z;
  k = domainnr(x, z, &alf, epiptr, gratptr);

  if (k == 0)
    {
      printf("\n (%lf %lf) is not in the interior domain\n", x, z);
      exit(1);
    }

  factor *= 1.0/alf;
  
  for (dmn=dmns, j=1; j<k; j++)
    {
      dmn = dmn->nextptr;
    }

  npnls = dmn->getnpanels();

  for (u=0.0, j=0; j<npnls; j++)
    {
      pnl = &(pnls[dmn->indp[j]]);
      pnlidxu = pnl->getidxu();
      pnlidxv = pnl->getidxv();
      i0 = pnlidxu;
      i1 = pnlidxv;
      kappa = dmn->getkap();

      calcp(pnl, xc, kappa, dmn->ornt[j], ptl);

      xcoll = pnl->getxcoll();
      nrm = pnl->getnrm();

      // this is to cancel the self-term if there is one
      inner = (x - xcoll[0])*nrm[0] + (z - xcoll[1])*nrm[1];

      if (fabs(inner) < EPS)
	ptl[1] = 0.0;

      type = pnl->gettype();

      if (type <= 1)
	{
	  u += double(dmn->ornt[j])*ptl[0]*sol0[i1] - ptl[1]*sol0[i0];
	}
      else if (type == 2)
	{
	  u -= lambda*( ptl[0]*sol0[i1] + ptl[1]*sol0[i0]);
	}
      else
	{
	  u += ptl[0]*sol0[i1] - ptl[1]*rhs[i1];
	}
    }
  return factor*u;
}
