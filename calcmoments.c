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
#include <cstdlib>

#include "structure.h"
#include "grating.h"
#include "gtoothdmn.h"

using namespace std;

extern "C" {
  // LU solve
  void  zgetrs_(char *trans, long *n, long *nrhs, complex<double> *a,
		long *lda, long *pvt, complex<double> *b,
		long *ldb, long *inf);
  // LU factorization
  void  zgetrf_(long *nrow, long *ncol, complex<double> *a, long *lda,
		long *pvt, long *inf);
}
int calcrhs(int order, int nrows, int npnls,
	    complex<double> *rhs, complex<double> *b,
	    grating *gratptr);

int calcrhs1(int order, int nrows, int npnls,
	    complex<double> *rhs, complex<double> *b,
	    grating *gratptr);

int calcp(gtpanel *p, double *x, complex<double> kappa, int ornt,
           complex<double> *ptl);

int calcfc(grating *gratptr, structure *epiptr, int nrows,
	   complex<double> *sol, complex<double> *mom);

int calcfc1(grating *gratptr, structure *epiptr, int nrows,
	   complex<double> *sol, complex<double> *mom);

int calcmoments(structure *epiptr, grating *gratptr, int nrows,
		complex<double> *lambda, complex<double> **at,
		complex<double> **sol0t, complex<double> **sol1t,
		complex<double> **rhst, complex<double> **a1t,
		complex<double> **bt, complex<double> ***mp,
		long **piv)
{
  int k;
  int ord2,ord4, order;
  int ncolsb, pnltype, npnls;
  int irow, ind, npts, sqord4;
  int i, j, p, ii;
  int dmnpnls, pnlidxu, pnlidxv;
  int posu, posv;
  long nrhs, lda, info;
  int nrowcolsb, nroword4;
  int icol, jj;

  double period;
  double *xc;

  long *ipiv;

  complex<double> gam;
  complex<double> kappa, orn;
  complex<double> ptl[CSIZE];
  complex<double> twoPi;
  complex<double> **mom;
  complex<double> *a1, *a, *b, *rhs, *sol0, *sol1;

  gtpanel *pnlptr; 
  gtpanel *pnls, *pnl;

  gtdomain *dmn;

  char trans;

  trans = 'N';
  twoPi = complex<double>(0.0,(2.0 * M_PI));
  ipiv = NULL;

  // get grating period
  period = gratptr->getperiod();

  // get beta guess
  gam = epiptr->getbetaguess();

  order = gratptr->getssspchrmncs();
  ord2 = (2 * order) + 1;
  ord4 = 2*ord2;

  if (gratptr->gzk == NULL)
    {
      gratptr->gzk = new complex<double>[ord2];
    }

  for (k=-order; k<=order; k++)
    {
      gratptr->gzk[k+order] = gam + ((static_cast<double>(k)*twoPi)/period);
    }

  // number of refined panels on type 3 and 4 surfaces == number of
  // cols in matrix B
  npnls = gratptr->gettotalpnls();
  pnlptr = gratptr->gtrefpnlptr;
  ncolsb = 0;

  while (pnlptr->nextptr != NULL)
    {
      pnltype = pnlptr->gettype();
      if (pnltype > 2)
	{
	  ncolsb++;
	}
      pnlptr = pnlptr->nextptr;
    }

  a = *at;
  b = *bt;
  a1 = *a1t;
  rhs = *rhst;
  sol0 = *sol0t;
  sol1 = *sol1t;
  
  ind = nrows * nrows;
  nrowcolsb = nrows * ncolsb;
  nroword4 = nrows * ord4;

  if (a == NULL)
    {
      a = new complex<double>[ind];
      a1 = new complex<double>[ind];
      b = new complex<double>[nrowcolsb];
      rhs = new complex<double>[nroword4];
      sol0 = new complex<double>[nroword4];
      sol1 = new complex<double>[nroword4];

      for (i=0; i<nroword4; i++)
	{
	  real(rhs[i]) = 0.0;
	  imag(rhs[i]) = 0.0;
	}
    }
  else
    {
      for (i=(ind-1); i>=0; i--)
	{
	  a[i] = 0.0;
	}
      for (i=(nrowcolsb-1); i>=0; i--)
	{
	  b[i] = 0.0;
	}
      for (i=(nroword4-1); i>=0; i--)
	{
	  real(rhs[i]) = 0.0;
	  imag(rhs[i]) = 0.0;
	}
    }

  if ((mom = *mp) == NULL)
    {
      mom = new complex<double>*[2];
      for (i=0; i<2; i++)
	{
	  sqord4 = ord4 * ord4;
	  mom[i] = new complex<double>[sqord4];
	}
      *mp = mom;
    }
  else
    {
      for (i=0; i<2; i++)
	{
	  for (ii=(ord4*ord4)-1; ii>=0; ii--)
	    {
	      mom[i][ii] = 0.0;
	    }
	}
    }

  // setup the BEM-stiffness matrix and its derivative wrt lambda
  // note that calcp() includes self-term
  for (irow = 0, dmn = gratptr->gtdmnptr; dmn != NULL; dmn = dmn->nextptr)
    {
      kappa = dmn->getkap();
      npts = dmn->getnpts();

      for (i = 0; i < npts; i++, irow++)
	{
	  pnls = gratptr->gtrefpnlptr;
	  for (p = 0; p < dmn->refindp[i]; p++)
	    {
	      pnls = pnls->nextptr;
	    }
	  xc = pnls->getxcoll();
	  dmnpnls = dmn->getnpanels();

	  for (j = 0; j < dmnpnls; j++)
	    {
	      pnl = gratptr->gtrefpnlptr;
	      for (p = 0; p < dmn->refindp[j]; p++)
		{
		  pnl = pnl->nextptr;
		}
	      calcp(pnl, xc, kappa, dmn->refornt[j], ptl);
	      pnlidxu = pnl->getidxu();
	      pnlidxv = pnl->getidxv();
	      posu = irow + (nrows * pnlidxu);
	      posv = irow + (nrows * pnlidxv);
	      pnltype = pnl->gettype();

	      if (pnltype == 2)
		{
		  ptl[0] *= *lambda;
		  ptl[1] *= *lambda;
		  a[posv] += ptl[0];
		  a[posu] += ptl[1];
		  period = gratptr->getperiod();
		  // derivative wrt gamma
		  a1[posv] = period * ptl[0];
		  a1[posu] = period * ptl[1];
		}
	      else
		{
		  orn = complex<double>(dmn->refornt[j],0);
		  a[posv] -= (orn * ptl[0]);  // ornt because of du/dn
		  a1[posv] = 0.0;
		  if (pnltype < 3)
		    {
		      a[posu] += ptl[1];
		      a1[posu] = 0.0;
		    }
		  else
		    {
		      b[posu] -= ptl[1];
		    }
		}

	    } // for panels
	} // for collocPts

    } // for domains

  calcrhs(order, nrows, npnls, rhs, b, gratptr);

  // LU factorization and solution
  if ((ipiv = *piv) == NULL)
    {
      ipiv = new long[nrows];
    }

  lda = static_cast<unsigned long> (nrows);
  nrhs = static_cast<unsigned long> (ord4);
  zgetrf_(&lda, &lda, a, &lda, ipiv, &info);

  if (info)
    {
       cout << "\n calcmoments(): zgetrf-info= " << info << endl;
       exit(1);
    }

  for (i=(nrows*ord4)-1; i>=0; i--)
    {
      sol0[i] = rhs[i];
    }
  zgetrs_(&trans, &lda, &nrhs, a, &lda, ipiv, sol0, &lda, &info);

  if (info)
    {
      cout << "\n calcmoments(): zgetrs-info= " << info << endl;
      exit(1);
    }

  calcfc(gratptr, epiptr, nrows, sol0, mom[0]);

  // Done with 0-th order moments, now calculate 1-st order moments.
  // By product rule do the deriv (wrt gamma) of Fourier coefficients,
  // then add up the deriv of the solution vector.

  calcfc1(gratptr, epiptr, nrows, sol0, mom[1]);

  // calculate deriv of rhs wrt gamma
  for (i=(nrows*ord4)-1; i>=0; i--)
    {
      rhs[i] = 0;
    }

  calcrhs1(order, nrows, npnls, rhs, b, gratptr);

  // subtract deriv of influence matrix times solution
  for (icol=k=0; k<ord4; icol+=nrows, k++)
    {
      for (i=0; i<nrows; i++)
	{
	  for (jj=j=0; j<nrows; j++, jj+=nrows)
	    {
	      rhs[icol+i] -= a1[i+jj]*sol0[icol+j];
	    }
	}
    }

  // solve linear system using prevous LU factors to get deriv of soln
  for (i=(nrows*ord4)-1; i>=0; i--)
    {
      sol1[i] = rhs[i];
    }
  zgetrs_(&trans, &lda, &nrhs, a, &lda, ipiv, sol1, &lda, &info); 

  if (info)
    {
      cout << " calcMoments(): zgetrs-info = " << info << endl;
      exit(1);
    }

  // add Fourier coefficients of deriv of soln to 1st-order moments
  calcfc(gratptr, epiptr, nrows, sol1, mom[1]);

  *at = a;
  *a1t = a1;
  *bt = b;
  *rhst = rhs;
  *sol0t = sol0;
  *sol1t = sol1;
  *piv = ipiv;

  return 0;
}
