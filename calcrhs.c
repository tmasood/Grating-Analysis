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

// calculate the RHS for 0-th order moments:
// the boundary condition on the left end first
// then the boundary condition on the right end

#include <iostream>
#include <complex>
#include "grating.h"
#include "gtoothpnl.h"

int calcrhs(int order, int nrows, int npnls,
	    complex<double> *rhs, complex<double> *b,
	    grating *gratptr)
{
  int k, ii, j, i, posu;
  int pnltype;
  int pnlidxu;

  complex<double> tmp;
  gtpanel *pnl;

  double *pnlxcoll;

  for (k=-order, ii=0; k<=order; k++, ii+=nrows)
    {
      pnl = gratptr->gtrefpnlptr;
      for (j=0; j<npnls; j++)
	{
	  pnltype = pnl->gettype();
	  
	  if (pnltype == 3)
	    {
	      pnlxcoll = pnl->getxcoll();
	      tmp = exp(gratptr->gzk[k+order] * pnlxcoll[1]);
	      pnlidxu = pnl->getidxu();
	      posu = nrows * pnlidxu;
	      for (i=0; i<nrows; i++)
		{
		  rhs[i+ii] += (b[i+posu]*tmp);
		}
	    }
	  
	  pnl = pnl->nextptr;
	}
    }

  for (k=-order; k<=order; k++, ii+=nrows)
    {
      pnl = gratptr->gtrefpnlptr;
      for (j=0; j<npnls; j++)
	{
	  pnltype = pnl->gettype();
	  
	  if (pnltype == 4)
	    {
	      pnlxcoll = pnl->getxcoll();
	      tmp = exp(gratptr->gzk[k+order] * pnlxcoll[1]);
	      pnlidxu = pnl->getidxu();
	      posu = nrows * pnlidxu;
	      for (i=0; i<nrows; i++)
		{
		  rhs[i+ii] += b[i+posu] * tmp;
		}
	    }
	  pnl = pnl->nextptr;
	}
    }
  return 0;
}
