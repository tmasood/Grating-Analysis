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

extern complex<double> calcsolint( gtpanel *pnls, gtdomain *dmns, double x,
				   double z, complex<double> lambda,
				   complex<double> *rhs,
				   complex<double> **sol0t,structure *epiptr,
				   grating *gratptr);

// Compute the solution in the interior domain and dump the results
// in two files, one for the real, the other for the imaginary part
//  Parameters:
//      order    order of harmonics in z
//      nPtsX    number of points where to print the solution X-direction
//      nPtsZ    number of points where to print the solution Z-direction
//      pnls     panels used for solution of interior problem
//      npnls    number of panels
//      NS       vector in nullspace of DtN operator (left and right side)

void printsolint(int order, int nptsx, int nptsz, int npnls, int nrows, 
		 gtdomain *dmns, gtpanel *pnls, complex<double> *ns,
		 complex<double> *rhs, complex<double> **sol0t,
		 complex<double> lambda, structure *epiptr,
		 grating *gratptr)
{
  int i, j;
  int ord2;
  int ord4;
  double x, z, x0int, x1int, dx, dz;
  double period;
  double *x0;
  complex<double> u;
  gtpanel *pnl;
  int pnltype;
  char fnam[50]; 
  ofstream fpre, fpim;

  strcpy(fnam, "solre");
  fpre.open(fnam);
  strcpy(fnam, "solim");
  fpim.open(fnam);

  ord2 = 2*order+1;
  ord4 = 2*ord2;

  calcrhsint(order, nrows, npnls, ns, rhs, sol0t, gratptr);
   
  // find the lower/upper bound of the interior interval
  pnls = gratptr->gtrefpnlptr;
  for (i=0; i<npnls; i++)
    {
      pnltype = pnl->gettype();
      x0 = pnl->getx0();

      if (pnltype == 3)
	{
	  x0int = x0[0];
	}
      if (pnltype == 4)
	{
	  x1int = x0[0];
	}
    }

  dx = (x1int - x0int)/nptsx;
  period = gratptr->getperiod();
  dz = period/nptsz;

  for (j=0, z=0.0; j<=nptsz; j++, z+=dz)
    {
      for (i=0, x=x0int; i<=nptsx; i++, x+=dx)
	{
	  u = calcsolint(pnls, dmns, x, z, lambda, rhs, sol0t, epiptr,
			 gratptr);
	  fpre << real(u) << " ";
	  fpim << imag(u) << " ";
	}
      fpre << endl;
      fpim << endl;
    }
  
  fpre.close();
  fpim.close();
}
