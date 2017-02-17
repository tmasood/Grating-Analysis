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


int laymat(double tr, dcomplex wx, dcomplex kc, dcomplex yn,
	   dcomplex *tm1, dcomplex *tm2, dcomplex *tm3)
{

  /* calculate field transformation matrix for any layer */
  /* --tr = normalized radian thickness, tr=k0*(x2 - x1) */
  /* --wx = Propagation coefficient in x, normalized to k0, */
  /* kx= k0*wx */
  /* kc= frequency factor, a multiplier of k0 */
  /* yn= transverse wave admit/imped factor, yx= wx*yn */
  /* tm1, tm2, tm3, matrix elements, tm4=tm1 */

  dcomplex zz, th, sn, qh, sc, tmp;

  /* phase thickness, cosine and sine */
  th = tr * wx * kc;
  *tm1 = zcos(th);
  sn = zsin(th);
  zz = wx * yn * sn;

  __real__ *tm2 = (-1.0) * __imag__ zz;
  __imag__ *tm2 = __real__ zz;

  if (zabsq(th) > 0.1)
    {
      zz = sn/(wx*yn);
    }

  else
    {
      /* for too small th use power series for sin(th)/th */
      qh = th * th;
      sc = A0 - qh * (A2 - qh * (A4 - qh *
				 (A6 - qh * (A8 - qh * A10))));
      tmp = (tr * kc * sc);
      zz = tmp/yn;
    }

  __real__ *tm3 = (-1.0) * __imag__ zz;
  __imag__ *tm3 = __real__ zz;

  return (0);

}
