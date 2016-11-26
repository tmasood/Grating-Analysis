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

// routines to calculate panel potentials

#include <iostream>
#include <complex>
#include "gtoothpnl.h"

extern "C" {
void zbesj_(double*, double*, double*, int*, int*, double[], 
            double[], int*, int*);
void zbesy_(double*, double*, double*, int*, int*, double[], 
            double[], int*, double[], double[], int*);
void zbesh_(double*, double*, double*, int*, int*, int*, double[], 
            double[], int*, int*);
}

/* Globals in this file */
extern double *xleg, *xlog, *wlog, *wleg; // Gauss abscissas/weights
extern int ordleg, ordlog;                // Gauss quadrature orders

extern double EPS;

// calculate potentials due a panel
//
// returns the vector ptl that contains the values
//   ptl[0]  single layer potential
//   ptl[1]  double layer potential, includes self-term for interior problem
//
//  parameters
//     p         panel
//     x         field point
//     ornt      orientation of panel normal:  -1 interior
//

int calcp(gtpanel *p, double *x, complex<double> kappa, int ornt,
	   complex<double> *ptl)
{
  complex<double> fval[2];
  int i;
  int morder = 2, kode = 1, itype = 1, ierr, nz;
  double fnu = 0.0, chr[2], chi[2], wrkr[2], wrki[2], zr, zi;
  double len, t, rr, frac, *tn, *nm, x0[2], r[2], *xcol;
  double *px0;
  double tan[2], nrm[2];

  // initialize stuff
  len = p->getlength();
  tn = p->gettan();
  nm = p->getnrm();
  xcol = p->getxcoll();

  for (i = 0; i < 2; i++)
    {
      ptl[i] = 0.0;
      tan[i] = tn[i];
      nrm[i] = ornt*nm[i];
      x0[i] = x[i] - xcol[i];
    }

  // square of the distance measured from node i
  rr  = (pow(x0[0],2.0)) + (pow(x0[1],2.0));

  if (rr < EPS)
    {
      morder = 1;
      // self terms: Gauss-Legendre for real part,
      // Gauss-quadrature with logarithmic weight for imag part
      // integral is only over half length

      for (i = 0; i < ordleg; i++)
        {
	  // Gauss-Legendre for Bessel-J
          t = xleg[i]*len*0.5;
          zr = real(kappa)*t;
          zi = imag(kappa)*t;
          zbesj_(&zr, &zi, &fnu, &kode, &morder, chr, chi, &nz, &ierr);
          fval[0] = complex<double>(chr[0], chi[0]);
          ptl[0] += wleg[i]*fval[0];
        }

      for (i = 0; i < ordlog; i++)
        {
	  // logarithmic quad rule for Bessel-Y
          t = xlog[i]*len*0.5;
          zr = real(kappa)*t;
          zi = imag(kappa)*t;
          zbesy_(&zr, &zi, &fnu, &kode, &morder, chr, chi, &nz,
		 wrkr, wrki, &ierr);
          fval[0] = complex<double>(-chi[0],chr[0]);
          ptl[0] += wlog[i]*fval[0];
        }
      
      ptl[0] *= len;            // normalize to interval length
      ptl[1] = complex<double>(0.0,-2.0);   // self-term for interior problem
    }
  else
    {
      // distant interactions: Gauss-Legendre quadrature
      px0 = p->getx0();
      x0[0] = x[0] - px0[0];
      x0[1] = x[1] - px0[1];
      // printf("x0[0] = %lf  x0[1] = %lf \n",x0[0], x0[1]);

      morder = 2;
      for (i = 0; i < ordleg; i++)
        {
          t = xleg[i] * len;
          r[0] = x0[0] - t*tan[0];
          r[1] = x0[1] - t*tan[1];

          // distance measured from node i
          rr = sqrt(pow(r[0],2.0) + pow(r[1],2.0));

          // kappa * rr = argument of the hankel function
          zr = real(kappa)*rr;
          zi = imag(kappa)*rr;

          // Use zero order Hankel function
          // hankel function order fnu = 0.0
          // chr = real part of result vector
          // chi = imag part of the result vector
          // morder = superscript of Hankel function = 2
          // itype = number of terms in the Hankel sequence
          // zr = real part of argument
          // zi = imag part of argument
          zbesh_(&zr, &zi, &fnu, &kode, &itype, &morder, chr, chi, &nz, &ierr);

          // fval = function value
	  fval[0] = complex<double>(chr[0],chi[0]);
          fval[1] = complex<double>(chr[1],chi[1]);

          // single layer
          // wleg = Jacobian
          // fval = fundamental solution
          // ptl = hankel function * jacobian
          ptl[0] += wleg[i] * fval[0];
          // double layer
          frac = (r[0]*nrm[0] + r[1]*nrm[1])/rr;
	  ptl[1] += wleg[i] * fval[1] * frac;
	}

      // add wavenumber to double layer potentials
      ptl[1] *= kappa;

      // normalize to interval length
      for (i = 0; i < 2; i++)
	{
	  ptl[i] *= len;
	}
    }

  return 0;
}
