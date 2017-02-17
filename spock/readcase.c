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
#include "case.h"

int readcase(char *strbuffer, struct CASE **caseptr)
{
  int day, month, year;
  char *caseflagptr;

  (*caseptr)->peropt = 1;

  caseflagptr = (char *)strstr(strbuffer, "KASE="); 
  if (caseflagptr != NULL)
    {
      sscanf(caseflagptr + 5, "%4d%2d%2d", &year, &month, &day);
      printf("File created: %d-%d-%d \n\n",day, month, year);
    }

  caseflagptr = (char *)strstr(strbuffer, "EPS1="); 
  if (caseflagptr != NULL)
    {
      sscanf(caseflagptr + 5, "%f", &(*caseptr)->eps1);
      printf("EPS1: Epsilon for convergence test ");
      printf("on root search: %g \n",(*caseptr)->eps1);
    }
	  
  caseflagptr = (char *)strstr(strbuffer, "EPS2="); 
  if (caseflagptr != NULL)
    {
      sscanf(caseflagptr + 5, "%f", &(*caseptr)->eps2);
      printf("EPS2: Epsilon for convergence test ");
      printf("on magnitude of eigen function: %g \n",(*caseptr)->eps2);
    }

  caseflagptr = (char *)strstr(strbuffer, "GAMEPS="); 
  if (caseflagptr != NULL)
    {
      sscanf(caseflagptr + 7, "%f", &(*caseptr)->gameps);
      printf("GAMEPS: Precision of gamma calculation %g \n",
	     (*caseptr)->gameps);
    }

  caseflagptr = (char *)strstr(strbuffer, "QZMR="); 
  if (caseflagptr != NULL)
    {
      sscanf(caseflagptr + 5, "%f", &(*caseptr)->qzmr);
      printf("QZMR: initial guess of real part of root %g \n",
	     (*caseptr)->qzmr);
      (*caseptr)->intqzmr = (*caseptr)->qzmr;
    }

  caseflagptr = (char *)strstr(strbuffer, "QZMI=");
  if (caseflagptr != NULL)
    {
      sscanf(caseflagptr + 5, "%f", &(*caseptr)->qzmi);
      printf("QZMI: initial guess of imaginary part of root %g \n",
	     (*caseptr)->qzmi);
      (*caseptr)->intqzmi = (*caseptr)->qzmi;
    }

  caseflagptr = (char *)strstr(strbuffer, "PRINTF=");
  if (caseflagptr != NULL)
    {
      sscanf(caseflagptr + 7, "%f", &(*caseptr)->printf);
      printf("PRINTF: output near and far-field data %g \n",
	     (*caseptr)->printf);
    }

  caseflagptr = (char *)strstr(strbuffer, "INITGS=");
  if (caseflagptr != NULL)
    {
      sscanf(caseflagptr + 7, "%f", &(*caseptr)->initgs);
      printf("INITGS: automatically generate input guess: %g \n",
	     (*caseptr)->initgs);
    }

  caseflagptr = (char *)strstr(strbuffer, "AUTOQW=");
  if (caseflagptr != NULL)
    {
      sscanf(caseflagptr + 7, "%f", &(*caseptr)->autoqw);
      printf("AUTOQW: GaAs automatic quantum well");
      printf(" width generation %g \n", (*caseptr)->autoqw);
    }

  caseflagptr = (char *)strstr(strbuffer, "NFPLT=");
  if (caseflagptr != NULL)
    {
      sscanf(caseflagptr + 6, "%f", &(*caseptr)->nfplt);
      printf("NFPLT: calculate near field: %g \n",
	     (*caseptr)->nfplt);
    }

  caseflagptr = (char *)strstr(strbuffer, "FFPLT=");
  if (caseflagptr != NULL)
    {
      sscanf(caseflagptr + 6, "%f", &(*caseptr)->ffplt);
      printf("FFPLT: calculate far field: %g \n",
	     (*caseptr)->ffplt);
    }

  caseflagptr = (char *)strstr(strbuffer, "DXIN=");
  if (caseflagptr != NULL)
    {
      sscanf(caseflagptr + 5, "%f", &(*caseptr)->dxin);
      printf("DXIN: calculate near field plot step: %g \n",
	     (*caseptr)->ffplt);
    }

  return 0;
}
