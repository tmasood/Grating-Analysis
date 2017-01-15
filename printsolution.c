//
//  Copyright 2015, 2016 Quantum Designs LLC, Taha Masood,
//  Johannes Tausch  and Jerome Butler
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

extern complex<double> calcsolnpt(int order, complex<double> *v,
				  complex<double> *w, complex<double> *gkb,
				  double h);
extern int translatevw(int order, complex<double> *v, complex<double> *w,
		       complex<double> *gkb, double width);

//  print the solution for z=0 into a file
//  Parameters:
//      order    order of harmonics in z
//      x0,x1    lower, upper bound of interval to print
//      npts     number of points where to print the solution
//      pnls     panels used for solution of interior problem
//      npnls    number of panels
//      nrows    number of unknowns in interior problem
//      vleft    vector in nullspace of DtN operator (left side)
//      vright   vector in nullspace of DtN operator (right side)

int printsolution(int order, double x0, double x1, int npts, int npnls, 
		  int nrows, gtpanel *pnls, complex<double> *vleft,
		  complex<double> *vright, complex<double> *extd2nleft,
		  complex<double> *extd2nright, grating *gratptr, 
		  structure *epiptr, int *sheetsleft, int *sheetsright,
		  complex<double> *sol0)
{
  int i, ii, j, k, nptsleft;
  int ord2;
  int ord4;
  int pnltype, pnlidxu;
  int ntrnsleft, ntrnsright;
  double x0int, x1int, x, dx, h, hupper, *xp;
  double *xx0;
  double period;
  double wleft, wright;
  double *pnlxcoll;
  double omega, wvl;
  double indx_r, indx_i;
  double loss;

  gtpanel *pnl;
  complex<double> *u, *v, *w, *gkb, tmp;
  complex<double> twopi;
  complex<double> gam;
  complex<double> I;
  complex<double> kapleft, kapright;
  complex<double> kapleft0, kapright0;

  layer *layleftptr, *layrightptr;

  char fnam[50]; 

  ofstream outfile;

  ord2 = (2*order) + 1;
  ord4 = 2*ord2;

  I = complex<double>(0.0,1.0);

  strcpy(fnam, "sol");
  outfile.open(fnam);

  wvl = epiptr->getwavelength();
  omega = (2 * M_PI)/wvl;

  twopi = complex<double>(0.0,(2.0 * M_PI));
  period = gratptr->getperiod();

  // get grating layer
  ntrnsleft = gratptr->getgl();
  ntrnsright = (epiptr->gettotallayers() - ntrnsleft - 2);

  // ntrnsleft is one less than the grating layer
  ntrnsleft--;  

  layleftptr = epiptr->layerptr;

  indx_r = layleftptr->getlayerindex();
  loss = layleftptr->getlayerloss();
  // calculate imaginary part of index from loss
  indx_i = ((loss * wvl)/(4 * M_PI * 10000));
  kapleft = complex<double>(indx_r, indx_i);
  kapleft *= omega;

  // infinite layer left
  kapleft0 = kapleft;

  for (i = 0; i < ntrnsleft; i++)
    {
      layleftptr = layleftptr->nextptr;
    }

  layrightptr = epiptr->layerptr;  
  while (layrightptr->nextptr != NULL)
    {
      layrightptr = layrightptr->nextptr;
    }

  indx_r = layrightptr->getlayerindex();
  loss = layrightptr->getlayerloss();
  // calculate imaginary part of index from loss
  indx_i = ((loss * wvl)/(4 * M_PI * 10000));
  kapright = complex<double>(indx_r, indx_i);
  kapright *= omega;

  // infinite layer right
  kapright0 *= kapright;

  pnl = gratptr->gtrefpnlptr;
  // find the lower/upper bound of the interior interval
  for (i = 0; i < npnls; i++)
    {
      pnltype = pnl->gettype();
      // grating period interfacing the left layer structure
      if (pnltype == 3)
	{
	  xx0 = pnl->getx0();
	  x0int = xx0[0];
	}
      // grating period interfacing the right layer structure
      if (pnltype == 4)
	{
	  xx0 = pnl->getx0();
	  x1int = xx0[0];
	}
      pnl = pnl->nextptr;
    }

  v = new complex<double>[ord2];
  w = new complex<double>[ord2];
  gkb = new complex<double>[ord2];

  dx = (x1 - x0)/npts;

  // calculate the solution on the left of the grating layer
  nptsleft = static_cast<int>(((x0int - x0)/dx) + 1);
  gam = epiptr->getbetaguess();

  if (nptsleft > 1)
    {
      u = new complex<double>[(nptsleft + ntrnsleft)];
      xp = new double[(nptsleft + ntrnsleft)];
      
      for (k = 0; k < ord2; k++)
	{
	  v[k] = vleft[k];
	  w[k] = -extd2nleft[k]*vleft[k];
	}
    
      for (x = x0int, ii = 0, j = ntrnsleft; j > 0; j--)
	{
	  for (k = -order; k <= order; k++)
	    {
	      indx_r = layleftptr->getlayerindex();
	      loss = layleftptr->getlayerloss();
	      // calculate imaginary part of index from loss
	      indx_i = ((loss * wvl)/(4 * M_PI * 10000));
	      kapleft = complex<double>(indx_r, indx_i);
	      kapleft *= omega;
	      gkb[k+order] = sqrt(pow((twopi*static_cast<double>(k)/period
				       + gam),2.0) + pow(kapleft,2.0));
	    }

	  wleft = layleftptr->getlayerthickness();

	  for (h = 0; h < wleft; h += dx, x -= dx, ii++)
	    {
	      u[ii] = calcsolnpt(order, v, w, gkb, -h);
	      xp[ii] = x;
	    }

	  translatevw(order, v, w, gkb, -wleft);
	  layleftptr = layleftptr->prevptr;
	}

      for (k = -order; k <= order; k++)
	{
	  gkb[k+order] = sqrt(pow((twopi*static_cast<double>(k)/period
				   + gam),2.0) + pow(kapleft0,2.0));
	  gkb[k+order] *= I*(static_cast<double>(sheetsleft[order+k]));

	}

      hupper = x - x0;

      for (h = 0; h < hupper; h += dx, x -= dx, ii++)
	{
	  for (tmp = 0.0, k = 0; k < ord2; k++)
	    {
	      tmp += v[k]*exp(gkb[k]*h);
	    }
	  
	  u[ii] = tmp;
	  xp[ii] = x;
	}

      for (ii--; ii >= 0; ii--)
	{
	  outfile <<  xp[ii] << "  " << real(u[ii]/u[0]);
	  outfile << "  " << imag(u[ii]/u[0]) << endl;
	}
    }
  // ======================================================

  // ======================================================
  // calculate the solution in the grating layer
  pnl = gratptr->gtrefpnlptr;
  for (i = 0; i < npnls; i++)
    {
      pnltype = pnl->gettype();
      pnlidxu = pnl->getidxu();

      // panel at the start of the grating period
      if (pnltype == 0)
	{
	  for (tmp = 0.0, k = 0, ii = 0; k < ord2; k++, ii += nrows)
	    {
	      tmp += sol0[ii+pnlidxu]*vleft[k];
	    }
	  for (k = 0; k < ord2; k++, ii += nrows)
	    {
	      tmp += sol0[ii+pnlidxu]*vright[k];
	    }
	  pnlxcoll = pnl->getxcoll();
	  outfile << pnlxcoll[0] << "  " <<  real(tmp/u[0]);
	  outfile << "  " <<  imag(tmp/u[0]) << endl;
	}
      pnl = pnl->nextptr;
    }
  // =======================================================

  // =======================================================
  // calculate the solution on the right of the grating layer
  if (x1int < x1)
    {
      for (k = 0; k < ord2; k++)
	{
	  v[k] = vright[k];
	  w[k] = extd2nright[k]*vright[k];
	}

      for (x = x1int, j = ntrnsright; j > 0; j--)
	{
	  indx_r = layrightptr->getlayerindex();
	  loss = layrightptr->getlayerloss();
	  // calculate imaginary part of index from loss
	  indx_i = ((loss * wvl)/(4 * M_PI * 10000));
	  kapright = complex<double>(indx_r, indx_i);
	  kapright *= omega;	  

	  for (k = -order; k <= order; k++)
	    {
	      gkb[k+order] = sqrt(pow((twopi*static_cast<double>(k)/period
				       + gam),2.0) + pow(kapright,2.0));
	    }

	  wright = layrightptr->getlayerthickness();
	  for (h = 0; h < wright; h+=dx, x+=dx)
	    {
	      // print calcsoln2;
	      tmp = calcsolnpt(order, v, w, gkb, h);
	      outfile << x << "  " << real(tmp/u[0]);
	      outfile << "  " <<  imag(tmp/u[0]) << endl;
	    }

	  translatevw(order, v, w, gkb, wright);
	  layrightptr = layrightptr->prevptr;
	}

      for (k = -order; k <= order; k++)
	{
	  gkb[k+order] = sqrt(pow((twopi*static_cast<double>(k)/period
				   + gam),2.0) + pow(kapright0,2.0));
	  gkb[k+order] *= I*static_cast<double>(sheetsright[order+k]); 
	}

      hupper = x1 - x;

      for (h = 0; h < hupper; h += dx, x += dx)
	{
	  for (tmp=0.0, k=0; k<ord2; k++)
	    {
	      tmp += v[k]*exp(gkb[k]*h);
	    }
	  outfile << x << "  " << real(tmp/u[0]);
	  outfile << "  " << imag(tmp/u[0]) << endl;
	}
    }

  outfile.close();

  return 0;
 
}

