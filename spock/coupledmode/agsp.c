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

dcomplex agsp(double x, double dmid, double g, int lup, int m,
	      dcomplex aup, dcomplex hup, struct STRUCT *structptr,
	      dcomplex bup)
{
  double p;
  dcomplex ztmp, ztmp1;
  dcomplex ci;
  
  p = 2 * (x - dmid)/g;

  if (p < 0.0)
    {
      p = 0.0;
    }
  else if (p > 1.0)
    {
      p = 1.0;
    }
  else
    {
      p = p;
    }

  __real__ ztmp =  acos(p);
  __imag__ ztmp =  0.0;
   
  __real__ ci = 0.0;
  __imag__ ci = 1.0;

  if (lup == 1)
    {
      ztmp1 = zsin(m*ztmp) * (aup*zexp(ci*hup*x)) * 
	(aup*zexp(ci*hup*x)); 
    }
  else if (lup == structptr->totallayers)
    {

      ztmp1 = zsin(m*ztmp) * (bup*zexp(-ci*hup*x)) * 
	(bup * zexp(-ci*hup*x)); 
    }
  else
    {
      ztmp1 = zsin(m*ztmp) * (aup*zexp(ci*hup*x)+bup*zexp(-ci*hup*x)) *
	(aup*zexp(ci*hup*x) + bup*zexp(-ci*hup*x)) * 
	(aup*zcos(hup*x) + bup*zsin(hup*x)) * (aup*zcos(hup*x) + 
					       bup*zsin(hup*x));
    }

  return (ztmp1);
}

