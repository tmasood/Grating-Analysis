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

int int4pt(double *x, double *f, int n, double *ans, double *er, int *fail)
{
  int ia, ib, i, j, k;
  double int1, e, s, c;
  double h1, h2, h3, h4, d1, d2, d3;
  double r1, r2, r3, r4;

  ia = 0;
  ib = n;

  if (ia == ib)
    {
      *ans = 0.0;
      *er = 0.0;
    }

  else
    {
      int1 = 0.0;
      e = 0.0;
      s = 0.0;
      c = 0.0;
      r4 = 0.0;
      if ((ia == (n-1)) && (n > 4))
	{
	  j = n - 2;
	}
      else if (ia > 2)
	{
	  j = ia;
	}
      else
	{
	  j = 3;
	}

      if ((ib == 1) && (n > 4))
	{
	  k = 3;
	}
      else if (n > (ib + 2))
	{
	  k = ib + 1;
	}
      else
	{
	  k = (n - 1);
	}

      for (i = j; i < k; i++)
	{
	  if (i == j)
	    {
	      h2 = x[j - 1] - x[j - 2];
	      if (h2 == 0.0)
		{
		  *fail = 3;
		  return (0);
		}
	      d3 = (f[j - 1] - f[j - 2]) / h2;
	      h3 = x[j] - x[j - 1];
	      if (h3 == 0.0)
		{
		  *fail = 3;
		  return (0);
		}
	      d1 = (f[j] - f[j - 1])/h3;
	      h1 = h2 + h3;
	      d2 = (d1 - d3)/h1;
	      h4 = x[j+1] - x[j];
	      if (h4 == 0.0)
		{
		  *fail = 3;
		  return (0);
		}
	      r1 = (f[j + 1] - f[j])/h4;
	      r2 = (r1 - d1)/(h4 + h3);
	      h1 = h1 + h4;
	      if (h1 == 0.0)
		{
		  *fail = 3;
		  return (0);
		}
	      r3 = (r2 - d2)/h1;
	      if (ia == 0.0)
		{
		  int1 = h2*(f[1] + h2*(d3/2.0 - h2*(d2/6.0 - 
						     (h2 + 2.0*h3)*r3/12.0)));
		  s = pow(-h2,3)*(h2*(3.0*h2 + 5.0*h4) + 10.0*h3*h1)/60.0;
		}
	    }
	  else
	    {
	      h4 = x[i+1] - x[i];
	      if (h4 == 0.0)
		{
		  *fail = 3;
		  return (0);
		}
	      r1 = (f[i + 1] -f [i])/h4;
	      r4 = h4 + h3;
	      if (r4 == 0.0)
		{
		  *fail = 3;
		  return (0);
		}
	      r2 = (r1 - d1)/r4;
	      r4 = r4 + h2;
	      if (r4 == 0.0)
		{
		  *fail = 3;
		  return (0);
		}
	      r4 = (r3 - d3)/r4;
	    }
	  if ((i <= ib) && (i > ia))
	    {
	      int1 = int1 + h3*((f[i] + f[i-1])/2.0 - h3*h3*(d2 + r2 +
			       (h2 - h4)*r3)/12.0);
	      c = pow(h3,3)*(2.0*h3*h3+5.0*(h3*(h4 + h2) + 2.0*h4*h2))/120.0;
	      e = e+(c+s)*r4;
	      if (i == j)
		{
		  s = 2*c + s;
		}
	      else
		{
		  s = c;
		}
	    }
	  else
	    {
	      e = e + r4*s;
	    }

	  if (i == k)
	    {
	      if (ib == n)
		{
		  int1 = int1 + h4*(f[n]-h4*((r1/2.0) + h4*(r2/6.0 +
			 ((2.0*h3) + h4)*r3/12.0)));
		  e = e - pow(h4,3)*r4*(h4*(3.0*h4 + 5.0*h2) +
					10.0 * h3*(h2+h3+h4))/60.0;
		}
	      if (ib >= (n - 1))
		{
		  e = e + s*r4;
		}
	    }
	  else
	    {
	      h1 = h2;
	      h2 = h3;
	      h3 = h4;
	      d1 = r1;
	      d2 = r2;
	      d3 = r3;
	    }
	}
      *ans = int1;
      *er = e;
    }
  return (0);
}
