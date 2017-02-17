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
 * This software is provided "as is" without expressed or implied warranty
 * to the extent permitted by applicable law.
]*/

#include <stdio.h>
#include <string.h>
#include "util.h"
#include "struct.h"
#include "layer.h"

extern int ingaasp_adachi(struct LAYER *, struct STRUCT *);
extern int ingaasp_henry(struct LAYER *, struct STRUCT *);
extern int algainas_adachi(struct LAYER *, struct STRUCT *);
extern int algainas_mondry(struct LAYER *, struct STRUCT *);

int calcindex(struct LAYER *headlayerptr, struct STRUCT *structptr)
{
  struct LAYER *tmplayerptr;

  tmplayerptr = headlayerptr;
  while (tmplayerptr != NULL)
    {
      if (tmplayerptr->matsys != 0)
	{
	  switch (tmplayerptr->matsys)
	    {
	    case 1:
	      ingaasp_adachi(tmplayerptr, structptr);
	      break;
	    case 12:
	      ingaasp_henry(tmplayerptr, structptr);
	      break;
	    case 3:
	      algainas_adachi(tmplayerptr, structptr);
	      break;
	    case 13:
	      algainas_mondry(tmplayerptr, structptr);
	      break;
	    default:
	      printf("tmplayerstart matsys = %d \n",tmplayerptr->matsys);
	      printf("\n MATSYS %d not found \n",tmplayerptr->matsys);
	      exit(1);
	      break;
	    }
	}
      tmplayerptr = tmplayerptr->nextptr;
    }
  return 0;
}
