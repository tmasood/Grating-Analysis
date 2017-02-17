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

extern int numzlp(double, int, double, int *, int *, struct CASE *);
extern int zincr(int, int, struct CASE *, struct STRUCT *);

int loopz(int m, struct STRUCT *structptr, struct CASE *caseptr)
{
  int n, i;
  int nerr, rdnpt, line;

  nerr = 0;

  numzlp(structptr->zfinv[m][1], structptr->ilz[m][1],
	 structptr->zinc[m][1], &n, &nerr, caseptr);
  if (nerr != 0)
    {
      return (0);
    }

  for (i = 0; i < n; i++)
    {
      /* increment loop variables, reset if first loop */
      zincr(i, m, caseptr, structptr);

      if (m >= structptr->loopzv)
	{
	  /* if outermost loop, then re-generate layer information   */
	  /* for outermost loop */
	  rdnpt = 0;
	  line = 1;
	  /* structs(line, rdnpt, by, tfile); 
	  ln = line;
	  if (nerr >= 0)
	    {
	      return (0);
	    }
	  */
	  /* generate initial guess (QZM) for current loop */
	  /* if (initgs > 0)
	    {
	      guess();
	    }
	  */
	  /* calculate single case */
	  /* 
	  j = 1;
	  loopx(i,j);
	  */
	}
      else
	{
	  /* jump to the next nested z loop */
	  /*
	  m = m + 1;
	  loopc(m);
	  */
	}
    }

  /* decrement loop counter to next outer loop */
  /*  m = m - 1; */
 
  return (0);
}
