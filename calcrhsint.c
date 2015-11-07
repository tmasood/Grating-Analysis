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

// Calculate the solution (=field and flux) at the interfaces using
// the solution from the last interior D2N map calculation. 
// On panels of type 0 and 1, the result is stored in the vector sol, i.e.
//   type=0 or 1 =>     u( pnl->xColl ) = Sol[ pnl->idxU ];
//                  du/dn( pnl->xColl ) = Sol[ pnl->idxV ];
// On the interface between exterior and 
// interior domains (i.e., panels of type 3 and 4).
// the result is stored in the vector Rhs, namely on index positions which
// belong to the FLUX position of the panel. That is, if 
//   type=3 or 4 =>   u( pnl->xColl ) = Rhs[ pnl->idxV ];
//   but          du/dn( pnl->xColl ) = Sol[ pnl->idxV ];
// The solution of type 2 panels can be obtained by multiplying the
// corresponding
// type 0 panel with the factor lambda, i.e.
// type=2  =>       u( pnl->xColl ) = lambda*Sol[ pnl->idxU ];
//                du/dn( pnl->xColl ) = lambda*Sol[ pnl->idxV ];
// and is therefore not calculated by this routine.

int calcrhsint(int order, int nrows, int npnls, complex<double> *ns,
		complex<double> *rhs, complex<double> **sol0t, 
		grating *gratptr)
{
  int k, i, i0, i1, icol;
  int ord2, ord4;
  int pnltype;

  double *xc;

  complex<double> tmp, u, v;
  complex<double> *sol0;
  gtpanel *pnl;

  ord2 = 2*order+1;
  ord4 = 2*ord2;

  sol0 = *sol0t;

  pnl = gratptr->gtrefpnlptr;
  for (i = 0; i < npnls; i++)
    {
      i0 = pnl->getidxu();
      i1 = pnl->getidxv();

      pnltype = pnl->gettype();
      if (pnltype <= 1)
	{
	  for (u=v=0.0, icol=0, k=0; k<ord4; icol+=nrows, k++)
	    {
	      u += ns[k]*sol0[icol + i0];
	      v += ns[k]*sol0[icol + i1];
	    }
	  sol0[i0] = u;
	  sol0[i1] = v;
	}

      else if (pnltype == 3)
	{
	  for (v=0.0, icol=0, k=0; k<ord4; icol+=nrows, k++)
	    {
	      v += ns[k] * sol0[icol+i1];
	    }

	  sol0[i1] = v;
	  for (u=0.0, k=0; k<ord2; k++)
	    {
	      xc = pnl->getxcoll();
	      tmp = exp(gratptr->gzk[k]*xc[1]);
	      u += ns[k]*tmp;
	    }
	  rhs[i1] = u;
	}

      else if (pnltype == 4)
	{
	  for (v=0.0, icol=0, k=0; k<ord4; icol+=nrows, k++ )
	    {
	      v += ns[k] * sol0[icol+i1];
	    }
	  sol0[i1] = v;
	  
	  for (u = 0.0, k=0; k<ord2; k++)
	    {
	      xc = pnl->getxcoll();
	      tmp = exp(gratptr->gzk[k]*xc[1]);
	      u += ns[k+ord2]*tmp;
	    }
	  rhs[i1] = u;
	}
      pnl = pnl->nextptr;
    }

  *sol0t = sol0;
  return 0;
}
 
