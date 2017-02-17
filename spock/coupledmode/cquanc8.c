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

int cquanc8(int fun, double a, double b, double abserr, double relerr,
	    dcomplex result, double errest, int nofun, double flag,
	    double dmid, double g, int lup, int m, dcomplex aup,
	    dcomplex hup, struct STRUCT *structptr, dcomplex bup, 
	    int lbelow, dcomplex bbelow, dcomplex hbelow, dcomplex abelow)
{
  double w0,w1,w2,w3,w4,stone,temp;
  double esterr,tolerr,step,x[17],xsave[8][30];
  dcomplex qprev, qnow, qdiff, qleft, area, cor11;
  dcomplex qright[31], f[17], fsave[8][30];
  int levmin, levmax, levout, nomax, nofin, lev, nim, i, j;

  dcomplex CZERO;

  __real__ CZERO = 0.0;
  __imag__ CZERO = 0.0;

  /* ***   stage 1 ***   general initialization set constants. */

  levmin = 1;
  levmax = 30;
  levout = 6;
  nomax = 5000;
  nofin = nomax - 8*(levmax - levout + pow(2,(levout + 1)));

  /* trouble when nofun reaches nofin */

  w0 =   3956.0 / 14175.0;
  w1 =  23552.0 / 14175.0;
  w2 =  -3712.0 / 14175.0;
  w3 =  41984.0 / 14175.0;
  w4 = -18160.0 / 14175.0;

  /* initialize running sums to zero */

  flag = 0.0;
  result = CZERO;
  cor11  = CZERO;
  errest = 0.0;
  area   = CZERO;
  nofun = 0;
  if (a == b)
    {
      return (0);
    }

  /* ***   stage 2 ***   initialization for first interval */

  lev = 0;
  nim = 1;
  x[0] = a;
  x[16] = b;
  qprev  = CZERO;
  if (fun = 0)
    {
      f[0] = agsp(x[0], dmid, g, lup, m, aup, hup, structptr, bup);
    }
  else
    {
      f[0] = bgsp(x[0], dmid, g, lbelow, structptr, m, bbelow, hbelow,
		  abelow);
    }

  stone = (b - a) / 16.0;
  x[8] = (x[0]  + x[16]) / 2.0;
  x[4] = (x[0]  + x[8]) / 2.0;
  x[12] = (x[8] + x[16]) / 2.0;
  x[2] = (x[0]  + x[4]) / 2.0;
  x[6] = (x[4] + x[8]) / 2.0;
  x[10] = (x[8]  + x[12]) / 2.0;
  x[14] =  (x[12] + x[16]) / 2.0;
  for (j = 2; j < 17; j = j + 2)   
    {
      if (fun == 0)
	{
	  f[j] = agsp(x[j], dmid, g, lup, m, aup, hup, structptr, bup);
	}
      else
	{
	  f[j] = bgsp(x[j], dmid, g, lbelow, structptr, m, bbelow, hbelow,
		  abelow);
	}
    }
  nofun = 9;

do
  {
    do
      {
	/* ***   stage 3 ***   central calculation */
	/* requires qprev,x0,x2,x4,...,x16,f0,f2,f4,...,f16. */
	/* calculates x1,x3,...x15, f1,f3,...f15,qleft, */
	/* qright,qnow,qdiff,area. */

	x[1] = (x[0] + x[2]) / 2.0;
	if (fun == 0)
	  {
	    f[1] = agsp(x[1], dmid, g, lup, m, aup, hup, structptr, bup);
	  }
	else
	  {
	    f[1] = bgsp(x[1], dmid, g, lbelow, structptr, m, bbelow, hbelow,
		  abelow);
	  }
	for (j = 3; j < 16; j = j + 2)
	  {
	    x[j] = (x[j-1] + x[j+1]) / 2.0;
	    if (fun == 0)
	      {
		f[j] = agsp(x[j], dmid, g, lup, m, aup, hup, structptr, bup);
	      }
	    else
	      {
		f[j] = bgsp(x[j], dmid, g, lbelow, structptr, m, bbelow,
			    hbelow, abelow);
	      }
	  }
     
	nofun = nofun + 8;
	step = (x[16] - x[0]) / 16.0;
	qleft = (w0*(f[0] + f[8]) + w1*(f[1] + f[7]) + w2*(f[2] + f[6]) +
		 w3*(f[3]+f[5]) + w4*f[4]) * step;
	qright[lev + 1] = (w0*(f[8] + f[16]) + w1*(f[9] + f[15]) + 
			   w2*(f[10] + f[14]) + w3*(f[11] + f[13]) +
			   w4*f[12]) * step;
	qnow = qleft + qright[lev + 1];
	qdiff = qnow - qprev;
	area = area + qdiff;

	/* ***   stage 4 *** interval convergence test */
	esterr = cdabs(qdiff) / 1023.0;
	tolerr = dmax1(abserr,relerr*cdabs(area)) * (step/stone);

	/****   stage 5   ***   no convergence */
	/* locate next interval. */
	nim = 2 * nim;
	lev = lev + 1;

	/* store right hand elements for future use. */
	for (i = 0; i < 8; i++)
	  {
	    fsave[i][lev] = f[i + 8];     ;
	    xsave[i][lev] = x[i + 8];
	  }

	/* assemble left hand elements for immediate use. */

	qprev = qleft;
	for (i = 0; i < 8; i++)
	  {
	    j = -i;
	    f[(2*j) + 18] = f[j + 9];
	    x[(2*j) + 18] = x[j + 9];
	  }
      } while (lev < levmin);
    
    if (lev >= levmax)
      {
	/* current level is levmax */
	flag = flag + 1.0;
      }
    else if (nofun > nofin)
      {
	/* ***   stage 6   ***   trouble section */
	/* number of function values is about to exceed limit. */
	
	nofin = 2 * nofin;
	levmax = levout;
	flag = flag + (b - x[0]) / (b - a);
      }
    else if (esterr <= tolerr)
      {
      }

    /* ***   stage 7   ***   interval converged */
    /* add contributions into running sums. */

    result = result + qnow;
    errest = errest + esterr;
    cor11  = cor11  + qdiff / 1023.0;

    /* locate next interval. */
    
    while (nim != 2*(nim/2))
      {
	nim = nim/2;
	lev = lev-1;
      }

    nim = nim + 1;
    if (lev > 0)
      {
	/* assemble elements required for the next interval. */
	qprev = qright[lev];
	x[0] = x[16];
	f[0] = f[16];
	for (i = 0; i < 8; i++)
	 {
	   f[2*i] = fsave[i][lev];     
	   x[2*i] = xsave[i][lev];
	 }
      }
  } while (lev > 0);

/*  *** stage 8 ***  finalize and return */
 result = result + cor11;

 /* make sure errest not less than roundoff level */

 if (errest == 0.0)
   {
     return (0);
   }

 temp = cdabs(result) + errest;

 while (temp == cdabs(result))
   {
     errest = 2.0*errest;
     temp = cdabs(result) + errest;
   }

  return (0);
}
