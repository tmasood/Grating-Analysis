/*[
 * Copyright 2004, 2005  Taha Masood
 *
 * Permission to use, copy, and distribute this software and its
 * documentation for any purpose with or without fee is hereby granted,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.
 *
 * Permission to modify the software is granted, but not the right to
 * distribute the complete modified source code.  Modifications are to
 * be distributed as patches to the released version.  Permission to
 * distribute binaries produced by compiling modified sources is granted,
 * provided you
 *   1. distribute the corresponding source modifications from the
 *    released version in the form of a patch file along with the binaries,
 *   2. add special version identification to distinguish your version
 *    in addition to the base release version number,
 *   3. provide your name and address as the primary contact for the
 *    support of your modified version, and
 *   4. retain our contact information in regard to use of the base
 *    software.
 * Permission to distribute the released version of the source code along
* with corresponding source modifications in the form of a patch file is
 * granted with same provisions 2 through 4 for binary distributions.
 *
 * This software is provided "as is" without express or implied warranty
 * to the extent permitted by applicable law.
]*/

#include <stdio.h>
#include <math.h>
#include "util.h"
#include "case.h"
#include "layer.h"
#include "struct.h"
#include "modcon.h"
#include "output.h"

extern int int4pt(double *, double *, int, double *, double *, int *);
extern int pfail(int );

int nrmlzf(double *s, int n, double *xxft, dcomplex *fyxft, 
	   struct CASE *caseptr)
{
  /* this routine normalizes the function fy(x) stored in the array */
  /* fyxft by evaluating the integral a = integral of */
  /* fy(x)*conj(fy(x))*dx, where fy(x) is a complex function and */
  /* conj(fy(x)) is its complex conjugate */

  dcomplex t;
  int i;
  int ifail;
  double a;
  double er;

  /* function int4pt requires at least 4 points in order to evaluate */
  /* the integral */
  ifail = 0;
  int4pt(xxft, s, n, &a, &er, &ifail);
  if (ifail > 0)
    {
      pfail(ifail);
    }

  if (a == 0.0)
    {
      printf("WARNING .... THe integral of s(x) = fy(x)*conj(fy(x)) = 0 \n");
      printf("Cannot normalize fy(x). \n");
      return (0);
    }

  /* compute a normalizing factor such that the integral of the */
  /* normalizing factor fy(x)*conj(fy(x))*dx = pfield */
  a = sqrt(caseptr->pfield/a);

  for (i = 0; i < n; i++)
    {
      fyxft[i] = fyxft[i] * a;
      t = fyxft[i];
      s[i] = pow(__real__ t, 2) + pow(__imag__ t, 2);
    }

  /* evaluate integral of normalized fy(x) * conj(fy(x)) */
  /* check validity of normalization. Integral should be = p */
  
  ifail = 0;
  int4pt (xxft, s, n, &a, &er, &ifail);
  if (ifail > 0)
    {
      pfail(ifail);
    }
  
  return (0);
}
