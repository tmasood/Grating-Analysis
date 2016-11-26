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
#include "newton.h"
#include "grating.h"

using namespace std;

const complex<double> ZZERO = complex<double>(0,0);

int calcfc1(grating *gratptr, structure *epiptr, int nrows,
	   complex<double> *sol, complex<double> *mom)
{
  int jj, ii, icol, k;
  int order, pnltype;
  int ord2, pnlidxv;
  int ord4;

  complex<double> tmp;

  gtpanel *pnlptr;

  double fac;
  double pnlen, period;
  double *pnlxcoll;

  order = gratptr->getssspchrmncs();
  ord2 = (2 * order) + 1;
  ord4 = 2 * ord2;

  for (jj = icol=0; jj < (pow(ord4,2.0)); icol += nrows, jj += ord4)
    {
      pnlptr = gratptr->gtrefpnlptr;
      while (pnlptr->nextptr != NULL)
	{
	  pnltype = pnlptr->gettype();
	  if (pnltype > 2)
	    {
	      if (pnltype == 3)
		{
		  for (k=-order, ii=0; k<=order; k++, ii++)
		    {
		      pnlxcoll = pnlptr->getxcoll();
		      tmp = -pnlxcoll[1] * exp(-gratptr->gzk[k+order]
					       * pnlxcoll[1]);
		      pnlidxv = pnlptr->getidxv();
		      pnlen = pnlptr->getlength();
		      mom[ii+jj] += sol[icol + pnlidxv] * tmp * pnlen;
		    }
		}
	      else
		{
		  for (k=-order, ii=ord2; k<=order; k++, ii++)
		    {
		      pnlxcoll = pnlptr->getxcoll();
		      tmp = -pnlxcoll[1] * exp(-gratptr->gzk[k+order]
					       * pnlxcoll[1]);
		      pnlidxv = pnlptr->getidxv();
		      pnlen = pnlptr->getlength();
		      mom[ii+jj] += sol[icol + pnlidxv] * tmp * pnlen;
		    }
		}
	    }
	  pnlptr = pnlptr->nextptr;
	}
    }

  // get grating period
  period = gratptr->getperiod();
  fac = 1.0/period;
  for (jj=(ord4*ord4)-1; jj>=0; jj--)
    {
      mom[jj] *= fac;
    }

  return 0;
}

