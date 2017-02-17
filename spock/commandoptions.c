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
#include "struct.h"
#include "layer.h"

extern int genplotindexprofile(char *, struct LAYER *, struct STRUCT *);
extern int genindexprofile(char *, struct LAYER *);
extern int createlayerfile(char *, struct LAYER *);

int commandoptions(char *filename, struct LAYER *headlayerptr,
		   struct STRUCT *structptr)
{
  struct LAYER *tmplayerptr;
  int commandindex;

  printf("filename %s \n",filename);
  printf("\n \n \n");
  printf("Enter 1 --  create layer file\n");
  printf("Enter 2 --  generate index profile\n");
  printf("Enter 3 --  generate and plot index profile using gnuplot \n");
  printf("Enter 4 --  calculate overlap integral \n");
  printf("> ");
  scanf("%d", &commandindex);

  switch(commandindex)
    {
    case 1:
      createlayerfile(filename, headlayerptr);
      break;
    case 2:
      genindexprofile(filename, headlayerptr);
      break;
    case 3:
      genplotindexprofile(filename, headlayerptr, structptr);
      break;
    case 4:
      break;
    default:
      break;
      exit(1);
    }
  return 0;
}
