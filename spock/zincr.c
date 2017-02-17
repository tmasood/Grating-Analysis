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

int zincr(int i, int m, struct CASE *caseptr, struct STRUCT *structptr)
{
  /* Increments z loop variables. Resets variables to initial */
  /* values if start of loop. Increments otherwise */

  int ii;

  for (ii = 1; ii <= structptr->lzcnt[m]; ii++)
    {
      switch (structptr->ilz[m][ii])
	{
	case 1:
	  if (i == 0)
	    {
	      caseptr->qzmr = caseptr->intqzmr;
	    }
	  else
	    {
	      caseptr->qzmr = caseptr->qzmr + structptr->zinc[m][ii];
	    }
	  break;

	case 2:
	  if (i == 0)
	    {
	      caseptr->qzmi = caseptr->intqzmi;
	    }
	  else
	    {
	      caseptr->qzmi = caseptr->qzmi + structptr->zinc[m][ii];
	    }
	  break;
	  
	}
    }

  return (0);
}
