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
#include "util.h"
#include "layer.h"
#include "charmatrix.h"
#include "struct.h"
#include "modcon.h"
#include "local_complex.h"

/* This function calculates propagation variables, wave */
/* admittance/impedance, phase thickness and integrals, */
/* characteristic matrices and eigen-equation function */
/* for a structure. */

/* To set the global variables an initial call */
/* to this function with kdo0 = 0 is required */

int charmatrix(struct LAYER *headlayerptr, struct STRUCT *structptr,
	       struct MODCON *modconptr, int kdo0)
{
  dcomplex CZERO, CONE;
  dcomplex wc, ac, bc, dc, cc;
  dcomplex qx1, qx2, q0;

  int kdid;
  struct LAYER *tmplayerptr;

  __real__ CZERO = 0.0;
  __imag__ CZERO = 0.0;
  __real__ CONE = 1.0;
  __imag__ CONE = 0.0;

  kdid = kdo0;
  tmplayerptr = headlayerptr->nextptr;

  /* calc variables for finite layers */
  while (tmplayerptr->nextptr != NULL)
    {
      /* transverse propagation coeff for inner layers */
      wc = zsqrt(tmplayerptr->qn - structptr->qz);
      tmplayerptr->wx = wc;
      tmplayerptr->yx = wc * tmplayerptr->yn;
      tmplayerptr = tmplayerptr->nextptr;
    }

  if (kdid <= 1)
    {
      return (0);
    }

  /* calculate phase thickness for each layer */
  /* calculate phase integral with real (phr) */
  /* and imaginary (img) part */
  structptr->phr = CZERO;
  structptr->phi = CZERO;

  /* calculate phase thickness and integral for finite layers. */
  /* headlayerptr points to the first infinite layer */
  tmplayerptr = headlayerptr->nextptr;
  /* while loop exit condition-- when tmplayerptr points to */
  /* the last infinite layer */
  while (tmplayerptr->nextptr != NULL)
    {
      ac = tmplayerptr->wx * structptr->kc * tmplayerptr->tr;
      tmplayerptr->ph = ac;
      structptr->phr = structptr->phr + ac;
      if (__imag__ ac < 0.0)
	{
	  ac = -ac;
	}
      structptr->phi = structptr->phi + ac;
      tmplayerptr = tmplayerptr->nextptr;
    }

  if (kdid <= 2)
    {
      return (0);
    }

  /* calculate layer matrices and system characteristic matrix */
  Tmatrix.m0 = CONE;
  Tmatrix.m1 = CZERO;
  Tmatrix.m2 = CZERO;
  Tmatrix.m3 = CONE;
  tmplayerptr = headlayerptr->nextptr;
  /* matrix for each layer */
  while(tmplayerptr->nextptr != NULL)
    {
      ac = tmplayerptr->ph;
      tmplayerptr->cl0 = zcos(ac);
      bc = zsin(ac);
      __real__ tmplayerptr->cl1 = -1.0 * __imag__ (bc * tmplayerptr->yx);
      __imag__ tmplayerptr->cl1 = __real__ (bc * tmplayerptr->yx);
      if (zabsq(ac) < 0.1)
	{
	  /* for too small a phase thickness use pwr series for */
	  /* sinc function */
	  dc = ac * ac;
	  cc = A0 - dc * (A2 - dc * (A4 - dc * (A6 - dc *
						(A8 - dc * A10))));
	  q0 = cc / tmplayerptr->yn;
	  ac = tmplayerptr->tr * structptr->kc * q0;
	  __real__ tmplayerptr->cl2 = -1.0 *__imag__ ac;
	  __imag__ tmplayerptr->cl2 = __real__ ac;
	}
      else
	{
	  __real__ tmplayerptr->cl2 = 
	    -1.0 * __imag__ (bc / tmplayerptr->yx);
	  __imag__ tmplayerptr->cl2 = __real__ (bc / tmplayerptr->yx);
	}
      /* accumulate system matrix */
      /* multiply by each layer matrix, from the left */
      ac = (tmplayerptr->cl0 * Tmatrix.m0) + 
	(tmplayerptr->cl2 * Tmatrix.m1);
      Tmatrix.m1 = (tmplayerptr->cl1 * Tmatrix.m0) + 
	(tmplayerptr->cl0 * Tmatrix.m1);
      Tmatrix.m0 = ac;
      bc = (tmplayerptr->cl0 * Tmatrix.m2) + 
	(tmplayerptr->cl2 * Tmatrix.m3);
      Tmatrix.m3 = (tmplayerptr->cl1 * Tmatrix.m2) + 
	(tmplayerptr->cl0 * Tmatrix.m3);
      Tmatrix.m2 = bc;
      tmplayerptr = tmplayerptr->nextptr;
    } /* end while loop */

  /* calculate determinant and antideterminant */
  cc = Tmatrix.m0 * Tmatrix.m3;
  dc = Tmatrix.m1 * Tmatrix.m2;
  structptr->detc[0] = cc - dc;
  structptr->detc[1] = cc + dc;

  if (kdid <= 3)
    {
      return (0);
    }

  /* impose boundry conditions at outer boundries */
  /* boundry conditions at first layer boundry. */
  /* calculate transverse wave propagation coefficient wx. Positive */
  /* real part. */
  ac = CONE;
  qx1 = headlayerptr->qn - structptr->qz;
  wc = zsqrt(qx1);
  /* change sign if necessary for principal branch half plane */
  /* of wx. */
  if ( __real__ (wc * zcnjg(modconptr->vpb1)) < 0.0)
    {
      wc = -wc;
    }
  headlayerptr->wx = wc;
  /* transverse wave admittance/impedance. yx for first outer layer */
  bc = wc * headlayerptr->yn;
  headlayerptr->yx = bc;

  /* save tangential field values at the first (outer) boundary. */
  /* Fy = 1 and Fz = -Y1 */
  headlayerptr->fy = ac;
  headlayerptr->fz = -bc;

  /* partial matrix product, from the right (top) by boundary matrix */
  cc = Tmatrix.m0 * ac;
  dc = Tmatrix.m2 * bc;
  Smatrix.m0 = cc + dc;
  Smatrix.m2 = cc - dc;
  cc = Tmatrix.m1 * ac;
  dc = Tmatrix.m3 * bc;
  Smatrix.m1 = cc + dc;
  Smatrix.m3 = cc - dc;

  /* boundry condition at second boundry (last layer)*/
  /* Transverse propagation coefficient wx */
  tmplayerptr = headlayerptr;
  while(tmplayerptr->nextptr != NULL)
    {
      tmplayerptr = tmplayerptr->nextptr;
    }
  /* positive real part of transverse propagation coefficient */
  qx2 = tmplayerptr->qn - structptr->qz;
  wc = zsqrt(qx2);
  /* Change sign if necessary for principal branch half-plane */
  /* of wx. */
  if ( __real__ (wc * zcnjg(modconptr->vpb2)) < 0.0)
    {
      wc = -wc;
    }
  tmplayerptr->wx = wc;
  /* transverse wave admit/imped, yx, for outer layer material */
  ac = wc * tmplayerptr->yn;
  tmplayerptr->yx = ac;
  bc = CONE;
  /* save tangential field values at the second (outer) boundary. */
  /* Fy = 1 and Fz = +Yn */
  tmplayerptr->fy = bc;
  tmplayerptr->fz = ac;
	  
  /* complete the matrix product, from the left (bottom) by inverse */
  /* boundary matrix */
  cc = ac * Smatrix.m0;
  dc = bc * Smatrix.m1;
  Smatrix.m0 = cc + dc;
  Smatrix.m1 = cc - dc;
  cc = ac * Smatrix.m2;
  dc = bc * Smatrix.m3;
  Smatrix.m2 = cc + dc;
  Smatrix.m3 = cc - dc;
  structptr->dets = 
    (Smatrix.m0 * Smatrix.m3) - (Smatrix.m1 * Smatrix.m2);

  if (kdid <= 4)
    {
      return (0);
    }
  
/* calculate logitudinal propagation variables, wave */
  /* admittance and impedance, and quadrant indicators */
  /* for qz and qx planes */
  structptr->wz = zsqrt(structptr->qz);
  tmplayerptr = headlayerptr;
  while (tmplayerptr != NULL)
    {
      tmplayerptr->yz = structptr->wz * tmplayerptr->yn;
      tmplayerptr = tmplayerptr->nextptr;
    }

  return 0;
}
