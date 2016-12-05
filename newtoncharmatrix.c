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

using namespace std;

const complex<double> ZZERO = complex<double>(0,0);

complex<double> zinv(complex<double> a);

// computes the solution of linear equations (a * x = b) using LU
// factorization.
extern "C" void zgesvx_(char *fact, char *trans, int *n, int *nhrs,
			complex<double> *a, int *lda, complex<double> *af,
			int *ldaf, int *ipiv, char *equed, double *r,
			double *c, complex<double> *b, int *ldb,
			complex<double> *x, int *ldx,  double *rcond,
			double *ferr, double *berr, complex<double> *work,
			double *rwork, int *info);

// zgeev_ computes the eigenvalues for the matrix "a"
extern "C" void zgeev_(char *jobvl, char *jobvr, int *n, complex<double> *a,
		       int *lda, complex<double> *w, complex<double> *vl,
		       int *ldvl, complex<double> *vr, int *ldvr,
		       complex<double> *work, int *lwork, double *rwork,
		       int *info);

matnew newtoncharmatrix(complex<double> *zg, structure *epiptr)
{
  int totallayers;
  int dim, dim2;
  int i, j, ij;
  int aij, bij, cij, dij;
  int info, lwork;

  layer *layerptr, *headlayptr;

  double dt;

  double *r, *c;
  double *ferr;
  double *berr;
  double *rwork;
  double rcond;
  double rmax, rval;
  double nlt;

  complex<double> *zt, *ztp, ct, st, ht;
  complex<double> hz0, hzll, hp, hd;
  complex<double> z2g;
  complex<double> *af, *x;
  complex<double> *work;
  complex<double> *zeign, zmin;
  complex<double> *zvl, *zvr;
  complex<double> c_index, c2_index;
  complex<double> arg, ac;
  complex<double> sinzh, coszh;


  char fact;
  char trans;
  char equed;
  char jobvl;
  char jobvr;
  
  int *ipiv;

  matnew newtonmat;

  totallayers = epiptr->gettotallayers();
  dim = (totallayers - 1) << 1;
  dim2 = dim * dim;

  zt = new complex<double>[dim2];
  ztp = new complex<double>[dim2];
  af = new complex<double>[dim2];
  x =  new complex<double>[dim2];
  work =  new complex<double>[dim2];
  zeign =  new complex<double>[dim2];
  zvl =  new complex<double>[dim2];
  zvr =  new complex<double>[dim2];

  ipiv = new int[dim];
  r = new double[dim];
  c = new double[dim];
  ferr = new double[dim];
  berr = new double[dim];
  rwork = new double[2*dim];

  z2g = pow(*zg,2.0);

  layerptr = epiptr->layerptr;
  headlayptr = epiptr->layerptr;

  while (layerptr != NULL)
    {
      c2_index = layerptr->getlayerkappa();
      nlt = layerptr->getlayerkth();
      // hz
      if ((headlayptr == layerptr) || (layerptr->nextptr == NULL))
	{
	  // the first layer or the last layer hz
	  ac = sqrt(-1.0 * (z2g + c2_index));
	  layerptr->sethz(ac);
	}
      else
	{
	  ac = sqrt(1.0 * (z2g + c2_index));
	  layerptr->sethz(ac);
	  arg = nlt * ac;
	  coszh = cos(arg);
	  layerptr->setczh(coszh);
	  sinzh = sin(arg);
	  layerptr->setszh(sinzh);
	}
      layerptr = layerptr->nextptr;
    }

  for (i = 0; i < dim2; i++)
    {
      zt[i] = ZZERO;
      ztp[i] = ZZERO;
    }

  layerptr = epiptr->layerptr;
  while (layerptr->nextptr->nextptr != NULL)
    {
      layerptr = layerptr->nextptr;
    }

  for(j = 0; j  < (totallayers - 2); j++)
    {  /* Fill in the Transfer matrices. */
      i = 2*j;
      aij = i + ((i+2)*dim);
      bij = i + ((i+3)*dim);
      cij = (i+1) + ((i+2)*dim);
      dij = (i+1) + ((i+3)*dim);

      ct = layerptr->getczh();
      zt[aij] = ct;
      st = layerptr->getszh();
      ht = layerptr->gethz();
      zt[bij] = st/ht;
      zt[cij] = -1.0 * st * ht;
      zt[dij] = ct;

      // negative derivatives so that det(ZT^{-1}*ZTP - \lambda*I) = 0
      hp = *zg/ht;
      dt = layerptr->getlayerkth();
      hd = ht * dt;
      ztp[aij] = dt * (st * hp);
      ztp[bij] = ((st - (hd * ct)) / pow(ht,2.0)) * hp;
      ztp[cij] = (st + (hd * ct)) * hp;
      ztp[dij] = ztp[aij];

      layerptr = layerptr->prevptr;
    }

  for(i = 0; i < (dim - 4); i++)
    {
      // fill in the unit matrices
      j = i + 4;
      ij = i + j*dim;
      zt[ij] = complex<double>(-1.0,0.0);
    }

  // fill in the odd elements
  zt[(3*dim)-2] = complex<double>(-1.0,0.0);
  zt[(4*dim)-1] = complex<double>(-1.0,0.0);

  zt[dim-4] = complex<double>(-1.0,0.0);
  // first layer hz
  hz0 = epiptr->layerptr->gethz();
  zt[dim-3] = hz0;

  zt[(2*dim)-2] = complex<double>(1.0,0.0);

  layerptr = epiptr->layerptr;
  while (layerptr->nextptr != NULL)
    {
      layerptr = layerptr->nextptr;
    }

  // last layer hz
  hzll = layerptr->gethz();
  zt[(2*dim)-1] = hzll;

  ztp[dim-3] = *zg / hz0;           /*negative derivative */
  ztp[(2*dim)-1] =  *zg / hzll;

  for (i = 0; i < dim; i++)
    {
      for (j = 0; j < dim; j++)
	{
	  // cout << zt[dim * i + j] << "  ";
	}
      // cout << endl << endl;
    }

  trans = 'N';
  fact  = 'E';
  equed = 'B';
  
  zgesvx_(&fact, &trans, &dim, &dim, zt, &dim, af, &dim, ipiv,
	  &equed, r, c, ztp, &dim, x, &dim, &rcond, ferr, berr,
 	  work, rwork, &info);

  cout <<  *zg << "         "  <<  rcond << "         " << info << endl;
  newtonmat.rcond = rcond;

  jobvl = 'N';
  jobvr = 'N';
  lwork = 2 * dim;

  // zeign - contains the computed eigenvales
  zgeev_(&jobvl, &jobvr, &dim, x, &dim, zeign, zvl, &dim,
           zvr, &dim, work, &lwork, rwork, &info);

  rmax = 0.0;
  for (i = 0; i < dim; i++)
    {
      rval = abs(zeign[i]);
      if (rval > rmax)
	{
	  rmax = rval;
	  zmin = zinv(zeign[i]);
	}
    }

  newtonmat.z = zmin;

  delete [] zt;
  delete [] ztp;
  delete [] af;
  delete [] x;
  delete [] work;
  delete [] zeign;
  delete [] zvl;
  delete [] zvr;

  return newtonmat;
}

complex<double> zinv(complex<double> a)
{
   double bm;
   double v_r;
   double v_i;
   complex<double> v;

   bm = (real(a)*real(a))+(imag(a)*imag(a));
   v_r = real(a)/bm;
   v_i =  (-1.0) * imag(a)/bm;
   v = complex<double>(v_r, v_i);
   return (v);
}
