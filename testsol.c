//
//  Copyright 2015, 2016 Quantum Designs LLC, Taha Masood, Johannes Tausch
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
#include "gtoothdmn.h"

using namespace std;

extern int calcrhsint(int order, int nrows, int npnls,
		      complex<double> *ns, complex<double> *rhs,
		      complex<double> **sol0t, grating *gratptr);

extern complex<double> calcsolint(gtpanel *pnls, gtdomain *dmns, double x,
				  double z, complex<double> lambda,
				  complex<double> *rhs,
				  complex<double> **sol0t,structure *epiptr,
				  grating *gratptr);

void testsol(int order, int nptsx, int nptsz, int npnls, int nrows, 
		 gtdomain *dmns, gtpanel *pnls, complex<double> *ns,
		 complex<double> *rhs, complex<double> **sol0t,
		 complex<double> lambda, structure *epiptr,
		 grating *gratptr)
{
  int j;
  int pnltype;
  int pnlidxu, pnlidxv;

  double x, z;
  double *pnlxcoll;

  complex<double> u, usol;
  complex<double> *sol0;

  gtpanel *pnl;

  calcrhsint(order, nrows, npnls, ns, rhs, sol0t, gratptr);
   
  for (j = 0; j < npnls; j++)
    {
      pnl = &(pnls[j]);
      pnlxcoll = pnl->getxcoll();
      x = pnlxcoll[0];
      z = pnlxcoll[1];
      u = calcsolint(pnls, dmns, x, z, lambda, rhs, sol0t, epiptr,
		     gratptr);

      pnltype = pnl->gettype();
      if (pnltype <= 1)
	{
	  pnlidxu = pnl->getidxu();
	  sol0 = *sol0t;
	  usol = sol0[pnlidxu];
	}
      else if (pnltype == 2)
	{
	  pnlidxu = pnl->getidxu();
	  usol = lambda*sol0[pnlidxu];
	}
      else 
	{
	  pnlidxv = pnl->getidxv();
	  usol = rhs[pnlidxv];
	}

      printf("(%lf,%lf) type=%d ex=(%lf,%lf) gr=(%lf,%lf) err=%lg\n",
	     x,z,pnltype, real(usol), imag(usol),
	     real(u), imag(u), abs(usol-u));
    }
}
