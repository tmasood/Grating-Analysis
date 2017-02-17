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

extern int nrmlzf(double *, int, double *, dcomplex *, struct CASE *);
extern int width(int, double *, double *, double *);
extern int int4pt(double *, double *, int, double *, double *, int *);
extern int pfail(int );

int ftfeld(char *filename, double *xxft, dcomplex *fyxft, dcomplex *pzxft,
	   int icnt, struct STRUCT *structptr, struct OUTPUT *outputptr,
	   struct CASE *caseptr, struct LAYER *headlayerptr)
{
  FILE *fpnf;
  FILE *fpff;

  char nffile[FILESIZE];
  char fffile[FILESIZE];

  dcomplex t, beta;
  double s[MAXPOINTS];
  double ta, tb, nfint, phase, thick, pzero;
  double rtodeg, temp, xft[MAXPOINTS], x, g0, i0;
  double nbarff, nbar2, btakon, er, b1, b2;
  double kosth, koxsth, itheta, theta;
  double ffield[MAXPOINTS], ffphas[MAXPOINTS];
  double dthetr, ffmax, sinth;
  double sinkxs[MAXPOINTS], coskxs[MAXPOINTS];
  double a1, a2, a3, a4, maxnfint;
  double costh, srn2st, gtheta, thetaa[MAXPOINTS];

  int i, j, n, ifail, ith;
  int ntotal, nlines;

  struct LAYER *tmplayerptr;

  rtodeg = 180.0/M_PI;

  strcpy(nffile, filename);
  strcat(nffile, ".nf");

  strcpy(fffile, filename);
  strcat(fffile, ".ff");

  fpnf = fopen(nffile,"w");
  fpff = fopen(fffile,"w");

  if (fpnf == NULL)
    {
      printf("Cannot open near field file to write \n");
      exit(1);
    }

  if (fpff == NULL)
    {
      printf("Cannot open far field file to write \n");
      exit(1);
    }

  /* Array xxft is assumed to be monotonically increasing. */
  /* Search for and eliminate any points (x, fy(x)) that */
  /* have any duplicate values. The latest duplicate value */
  /* is used in case of duplication. */

  temp = 1.0e3;
  n = 0;

  for (i = 0; i < icnt; i++)
    {

      t = fyxft[i];
      s[i] = pow(__real__ t, 2) + pow(__imag__ t, 2);
    }

  /* normalize the field fy(x) (stored in fyxft[]) by evaluating */
  /* the expression a = integral of fy(x) * conj(fy(x)*dx) and */
  /* then multiplying each fy(x) by sqrt(pfield/a). */
  nrmlzf(s,i, xxft, fyxft, caseptr);
  nrmlzf(s,i, xxft, pzxft, caseptr);
  width(i,xxft,s,&structptr->fwhpn);

  n = i;
  /* print near field intensity and phase */
  if (outputptr->nfout > 0)
    {
      /* write near field intensity and phase to the near field file */
      fprintf(fpnf,"xxft \t nfint \t phase \t nfreal \t nfimag \n");
    }

  maxnfint = 0.0;
  for (i = 0; i < icnt; i++)
    {
      maxnfint = MAX(maxnfint,fabs(__real__ pzxft[i]));
    }

  for (i = 0; i < icnt; i++)
    {
      t = fyxft[i];
      ta = __real__ t;
      nfint = fabs(__real__ pzxft[i])/maxnfint;
      phase = atan2(__imag__ t, __real__ t) * rtodeg;
      if (outputptr->nfout > 0 )
	{
	  fprintf(fpnf,"%g \t %g \t %g \t %g \t %g \n",
		  xxft[i], nfint, phase, ta, tb);
	}
    }

  /* compute the farfield. Set the x values so that zero is located */
  /* at the center of the device physically */
  thick = 0.0;
  tmplayerptr = headlayerptr->nextptr;
  while (tmplayerptr->nextptr != NULL)
    {
      thick = thick + tmplayerptr->tl;
      tmplayerptr = tmplayerptr->nextptr; 
    }

  pzero = thick/2.0;

  for (j = 0; j < i; j++)
    {
      xft[j] = xxft[j] - pzero;
    }

  /* compute some constant quantities and then evaluate the integrals */
  /* b1 and b2 */
  beta = structptr->wz * structptr->k0;
  nbarff = structptr->wz;
  nbar2 = pow(nbarff,2);
  btakon = beta / structptr->k0;

  /* form a table of real(fy(x)) and then evaluate the integral b1 */
  for (i = 0; i < n; i++)
    {
      s[i] = __real__ fyxft[i];
    }

  ifail = 0;
  int4pt(xft, s, n, &b1, &er, &ifail);
  if (ifail > 0)
    {
      pfail(ifail);
    }

  /* form a table of imag(fy(x)) and then evaluate the integral b2 */
  for (i = 0; i < n; i++)
    {
      s[i] = __imag__ fyxft[i];
    }

  ifail = 0;
  int4pt(xft, s, n, &b2, &er, &ifail);
  if (ifail > 0)
    {
      pfail(ifail);
    }

  /* compute g0 and i0 for theta = 0 */
  g0 = 2.0 * (nbarff + btakon) / (1.0 + nbarff);

  i0 = ((g0 * g0) * ((b1 * b1) + (b2 * b2)));

  /* f(theta) = 1 for theta = 0 */
  theta = 0.0;
  ffield[caseptr->ntheta] = i0;
  ffphas[caseptr->ntheta] = atan2(b2, b1) * rtodeg;
  dthetr = caseptr->dtheta / rtodeg;
  ffmax = i0;

  /* start of loop over theta ( 0 < theta >= 90 deg) */
  for (ith = 0; ith < caseptr->ntheta; ith++)
    {
      theta = theta + dthetr;
      sinth = sin(theta);
      kosth = structptr->k0*sinth;

      /* generate complete tables of cos(k0 * x * sin(theta)) and */
      /* sin(k0 * x * sin(theta)) for a full range of x values */
      for (i = 0; i < n; i++)
	{
	  x = xft[i];
	  koxsth = kosth*x;
	  sinkxs[i] = sin(koxsth);
	  coskxs[i] = cos(koxsth);
	  /* form a table of the quantity real(fy(x) *
	     cos(k0 * x * sin(theta)) */
	  s[i] = (__real__ fyxft[i]) * coskxs[i];
	}
      
      /* Integrate over x to get the integral a1 */
      ifail = 0.0;
      int4pt(xft, s, n, &a1, &er, &ifail);
      if (ifail > 0)
	{
	  pfail(ifail);
	}
      for (i = 0; i < n; i++)
	{
	  s[i] = (__real__ fyxft[i]) * sinkxs[i];
	}

      /* integrate over x to get the integral a2 */
      ifail = 0;
      int4pt(xft, s, n, &a2, &er, &ifail);
      if (ifail > 0)
	{
	  pfail(ifail);
	}
      for (i = 0; i < n; i++)
	{
	  s[i] = (__imag__ fyxft[i]) * coskxs[i];
	}

      /* integrate over x to get the integral a3 */
      ifail = 0;
      int4pt(xft, s, n, &a3, &er, &ifail);
      if (ifail > 0)
	{
	  pfail(ifail);
	}
      for (i = 0; i < n; i++)
	{
	  s[i] = (__imag__ fyxft[i]) * sinkxs[i];
	}

      /* integrate over x to get the integral a4 */
      ifail = 0;
      int4pt(xft, s, n, &a4, &er, &ifail);
      if (ifail > 0)
	{
	  pfail(ifail);
	}

      /* compute g(theta) */
      /* g(theta) is the obliquity factor. See paper "Radiation from */
      /* a solid-state laser", G. A. Hockham, July 1973, pg 389 - 391 */
      costh = cos(theta);
      srn2st = sqrt((nbar2 - pow(sinth,2)));
      gtheta = 2.0 * costh * (srn2st + btakon) / (costh + srn2st);

      /* finally compute i(theta) and the far field f(theta) */
      itheta = pow(gtheta,2) * (pow((a1 - a4),2) + pow((a2+a3),2));

      /* normalize far-field so that the maximum value is 1 */
      ffield[caseptr->ntheta+ith+1] = itheta;
      ffmax = MAX(ffmax, itheta);
      ffphas[caseptr->ntheta+ith+1] = atan2((a2+a3),(a1-a4))*rtodeg;
  
      /* To get f(-theta), reverse the sign of each sin(theta) term but */
      /* keep the sign of each cos(theta) term since sin(-u) = -sin(u) */
      /* but cos(-u) = cos(u) */
      a2 = -a2;
      a4 = -a4;
      itheta = pow(gtheta,2) * (pow((a1 - a4),2) + pow((a2+a3),2));
      ffield[caseptr->ntheta-ith+1] = itheta;
      ffmax = MAX(ffmax, itheta);
      ffphas[caseptr->ntheta-ith+1] = atan2((a2+a3),(a1-a4))*rtodeg;
    }

  /* when the entire array ffield has been generated, it should */
  /* contain (2*ntheta+1) values of f(theta) for increasing */
  /* values of theta from -thetax to +thetax in increments of */
  /* dtheta. Thus ffield(i), (where i = 1, 2, 3, ........, ntheta) */
  /* = f(theta). For negative values of theta (from -thetax, -thetax */
  /* + dtheta, -thetax + 2*dtheta, .....-dtheta) ffield(ntheta+1) = */
  /* f(theta) and for theta = 0 ffield(ntheta+1+i). */

  /* write far field data */
  ntotal = (2 * caseptr->ntheta) + 1;
  if (outputptr->ffout > 0)
    {
      fprintf(fpff,"theta \t ffield \t ffphase \n");
    }

  theta = -caseptr->thetax;

  for (i = 0; i < ntotal; i++)
    {
      if (abs(theta) < (1.0e-4 * caseptr->dtheta))
	{
	  theta = 0;
	}
      ffield[i] = ffield[i]/ffmax;
      if(outputptr->ffout > 0)
	{
	  fprintf(fpff,"%g \t %g \t %g \n", theta, ffield[i], ffphas[i]);
	}
      thetaa[i] = theta;
      theta = theta + caseptr->dtheta;
    }

  /* calculate full width half power of the far field */
  width(ntotal, thetaa, ffield, &structptr->fwhpf);

  /* print table of x, near field, theta, and far field intensity */
  nlines = MAX(n, ntotal);

  for (i = 0; i < nlines; i++)
    {
      if (i > n)
	{
	  break;
	}
      /* list near field intensity and phase */
      nfint = pow((__real__ fyxft[i]),2) + pow((__imag__ fyxft[i]),2);
      t = fyxft[i];
      phase = 0.0;
      if (abs(t) != 0.0)
	{
	  phase = atan2(__imag__ t, __real__ t) * rtodeg;
	}
      theta = theta + caseptr->dtheta;
    }
  fclose(fpff);
  fclose(fpnf);
  return (0);
}
