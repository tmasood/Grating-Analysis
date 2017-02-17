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
#include <string.h>
#include "util.h"
#include "case.h"
#include "layer.h"
#include "struct.h"
#include "modcon.h"
#include "output.h"
#include <math.h>

int guideparams(struct CASE *caseptr, struct STRUCT *structptr,
	   struct LAYER *headlayerptr, struct MODCON *modconptr,
		struct OUTPUT *outputptr, char *filename)
{
  FILE *fpout;

  struct LAYER *tmplayerptr;
  double nimag;
  dcomplex CONE, CZERO;
  dcomplex betaguess;
  dcomplex dc, bc, cc, ac;

  double tm;

  int lm, ll, lk, layernum, i;

  char outfile[FILESIZE];
  
  lm = ll = lk = 0;

  __real__ CONE = 1.0;
  __imag__ CONE = 0.0;
  __real__ CZERO = 0.0;
  __imag__ CZERO = 0.0;

  strcpy(outfile, filename);
  strcat(outfile, ".out");
  
  printf(" \n \n \t \t SPOCK \n \n");
  printf(" For calculation of complex eigen-modes \n");
  printf(" in multilayered structures \n \n \n");

  /* ***************************************************** */
  __real__ structptr->wz = 0.0;
  __imag__ structptr->wz = 0.0;

  __real__ structptr->qz = 0.0;
  __imag__ structptr->qz = 0.0;

  __real__ betaguess = caseptr->qzmr;
  __imag__ betaguess = caseptr->qzmi;
  caseptr->qzm[0] = betaguess;

  caseptr->kdof = 2;
  caseptr->kouf = -1;
  caseptr->pfield = 1;
  outputptr->nfout = 1;
  outputptr->ffout = 1;

  for (i = 0; i < structptr->mn; i++)
    {
      caseptr->phm[i] = CZERO;
      tmplayerptr = headlayerptr;
      while (tmplayerptr != NULL)
	{
	  tmplayerptr->ra[i] = 0.0;
	  tmplayerptr = tmplayerptr->nextptr;
	}
    }

  /* ***************************************************** */
  if (caseptr->lxyopt >= 0)
    {
      if (caseptr->kase != caseptr->last)
	{
	  caseptr->ksub = 0;
	  caseptr->last = caseptr->kase;
	}
      caseptr->ksub++;
      if (outputptr->modout)
	{
	  fpout = fopen(outfile,"w");
	  if (fpout == NULL)
	    {
	      printf("Cannot open output file to write \n");
	      exit(1);
	    }
	  fprintf(fpout,"\n **SPOCK (Simulation tool for Photonics and ");
	  fprintf(fpout,"Optical Circuit Knowledge)** \n");
	  fprintf(fpout,"===================================================");
	  fprintf(fpout,"============ \n");
	  fprintf(fpout,"Version: 0.01 \n\n\n");
	  fprintf(fpout,"** case number %d  %d ** \n\n",
		  caseptr->kase, caseptr->ksub);
	}
      
      if (headlayerptr == NULL)
	{
	  printf("ERROR:-There must be atleast one layer in your");
	  printf("structure \n");
	  exit(1);
	}

      if (structptr->mn <= 0)
	{
	  printf("ERROR:-Program must search for atleast one mode \n");
	  exit(1);
	}
      if (caseptr->lxyopt != 0)
	{
	  if (structptr->ny <= 0)
	    {
	      printf("ERROR:- Number of y-slices must be greater than 0 \n");
	      exit(1);
	    }
	  if (structptr->mny <= 0)
	    {
	      printf("ERROR:- Program must search for atleast one y-mode \n");
	      exit(1);
	    }
	}
    }

  /* dielectric permittivity */
  if (caseptr->peropt == 1)
    {
      structptr->k0 = (2 * M_PI)/structptr->wvl;
    }

  cc = structptr->kf * structptr->kf;

  tmplayerptr = headlayerptr;
  while(tmplayerptr != NULL)
    {
      nimag = tmplayerptr->nloss/structptr->k0;
      tmplayerptr->per = pow(tmplayerptr->nreal,2) - pow(nimag,2);
      tmplayerptr->pei = 2.0 * tmplayerptr->nreal * nimag;
      tmplayerptr->pmr = 1.0;
      tmplayerptr->pmi = 0.0;
      __real__ dc = tmplayerptr->per;
      __imag__ dc = tmplayerptr->pei;
      __real__ bc = tmplayerptr->pmr;
      __imag__ bc = tmplayerptr->pmi;
      tmplayerptr->qn = cc * bc * dc;

      /* kpol = 0 is for TE mode */
      /* kpol = 1 for TM mode */
      if (modconptr->kpol == 0)
	{
	  tmplayerptr->yn = CONE / (structptr->kf * bc);
	}
      else if (modconptr->kpol == 1)
	{
	  tmplayerptr->yn = CONE / (structptr->kf * dc);
	}
      else
	{
	  printf("\n  kpol = 0 is for TE mode calculations, \n");
	  printf("\n  kpol = 1 is for TM mode calculations. \n");
	  printf("No other value for KPOL parameter is allowed.\n");
	  exit(1);
	}

      tmplayerptr->kappa = zsqrt(tmplayerptr->yn - betaguess);
      lm++;
      tmplayerptr = tmplayerptr->nextptr;
    }

  tmplayerptr = headlayerptr;
  if (tmplayerptr != NULL)
    {
      tmplayerptr->bl = 0.0;
      tmplayerptr->xl = 0.0;
      tmplayerptr->tr = 0.0;
    }

  while(tmplayerptr->nextptr != NULL)
    {
      tmplayerptr->nextptr->bl = tmplayerptr->bl + tmplayerptr->nextptr->tl;
      tmplayerptr->nextptr->xl = tmplayerptr->xl + tmplayerptr->nextptr->tl;
      tmplayerptr->nextptr->tr = structptr->k0 * tmplayerptr->nextptr->tl;
      tmplayerptr->nextptr->ph = tmplayerptr->nextptr->tr *
	tmplayerptr->nextptr->kappa;
      tmplayerptr = tmplayerptr->nextptr;
    }

  /* control parameter for eig, non-eigen condition. One */
  /* or both boundaries */
  modconptr->kbc0 = 0.0;
  if (modconptr->kbc1 >= 1)
    {
      modconptr->kbc0 = 1;
    }
  if (modconptr->kbc2 >= 1)
    {
      modconptr->kbc0 = modconptr->kbc0 + 2;
    }

  __real__ modconptr->vpb1 = zcos(modconptr->apb1 * M_PI);
  __imag__ modconptr->vpb1 = zsin(modconptr->apb1 * M_PI);
  __real__ modconptr->vpb2 = zcos(modconptr->apb2 * M_PI);
  __imag__ modconptr->vpb2 = zsin(modconptr->apb2 * M_PI);

  ll = lm - 1;
  lk = lm - 2;

  /* print summary of structure */
  if (outputptr->modout == 1)
    {
      cc = structptr->kc * structptr->k0;
      fprintf(fpout,"Summary of the structure, polarisations");
      fprintf(fpout," and boundary conditions \n");
      fprintf(fpout,"---------------------------------------");
      fprintf(fpout,"---------------------- \n \n");
      fprintf(fpout,"Implied unit of length - um (1.0e-6 meters). \n");
      fprintf(fpout,"All parameters are normalized. \n");
      fprintf(fpout,"Nominal wavelength = %f um. \n",structptr->wvl);
      fprintf(fpout,"Free space k0 = 2*PI/wvl =  %f rad/um. \n",structptr->k0);
      fprintf(fpout,"Complex freqency factor kc = (%f, %f). \n",structptr->kc);
      fprintf(fpout,"Effective-k0 = kc * k0 = (%f, %f). \n",cc);
      fprintf(fpout,"All propagation coefficients are normalized to ");
      fprintf(fpout,"effective-k0.\nRelative frequency kf = (%f, %f). \n\n",
	      structptr->kf);
      fprintf(fpout,"Polarization kpol (=0 TE case, =1 TM case) = %d. \n",
	      modconptr->kpol);
      if (modconptr->kpol == 0)
	{
	  fprintf(fpout,"Transverse electric case- Ex = Hy = Ez = 0.\n");
	  fprintf(fpout,"Tangential Fy = Ey, longitudinal Fz = Hz and \n");
	  fprintf(fpout,"transverse Fx = -Hx and Yx = Fz/Fy and \n");
	  fprintf(fpout,"Yz = Fx/Fy are wave admittances in the \n");
	  fprintf(fpout,"postive x and z direction \n");
	}
      else
	{
	  fprintf(fpout,"transverse magnetic case- Hx = Ey = Hz = 0.\n");
	  fprintf(fpout,"Tangential Fy = Hy, longitudinal Fz = -Ez and \n");
	  fprintf(fpout,"transverse Fx = Ex. Yx = Fz/Fy and \n");
	  fprintf(fpout,"Yz = Fx/Fy are wave admittances in the \n");
	  fprintf(fpout,"postive x and z direction \n");
	}
      fprintf(fpout," \n \n Layered Structure: \n");
      fprintf(fpout,"Total layers = %d. \n",lm);
      fprintf(fpout,"Boundries = %d. \t Finite layers = %d", ll, lk);
      fprintf(fpout,"  \n\n");

      tmplayerptr = headlayerptr;
      fprintf(fpout,"-----------------------------------------------\n\n");
      fprintf(fpout,"Semi-infinite Layer: \n");
      fprintf(fpout, "\n Outer boundary is boundry 1 \n");
      ac = zsqrt(tmplayerptr->qn);
      fprintf(fpout,"pe = (%f %f), pm = (%f %f) \n",
	      tmplayerptr->per, tmplayerptr->pei,
	      tmplayerptr->pmr, tmplayerptr->pmi);
      fprintf(fpout,"layer number = 1 \n");
      fprintf(fpout,"qn = (%f %f), yn = (%f %f) \n",
	      tmplayerptr->qn, tmplayerptr->yn);
      fprintf(fpout,"rn = (%f %f) \n\n", ac);
     fprintf(fpout,"-----------------------------------------------\n\n");

      tmplayerptr = headlayerptr->nextptr;
      layernum = 2;
      fprintf(fpout,"Finite Layers: \n");
      while (tmplayerptr->nextptr != NULL)
	{
	  tm = tmplayerptr->tl / structptr->wvl;
	  bc = zsqrt(tmplayerptr->qn);
	  cc = tmplayerptr->tr * bc * structptr->kc;

	  fprintf(fpout,"pe = (%f, %f), pm = (%f, %f). \n",
		  tmplayerptr->per, tmplayerptr->pei,
		  tmplayerptr->pmr, tmplayerptr->pmi);
	  fprintf(fpout,"layer number = %d, tl = %f. \n", layernum,
		  tmplayerptr->tl);
	  fprintf(fpout,"qn = (%f, %f), yn = (%f, %f). \n",
		  tmplayerptr->qn, tmplayerptr->yn);
	  fprintf(fpout,"tl/wvl = %f, rn = (%f %f), ph0 = (%f %f). \n",
		  tm, bc, cc);
	  fprintf(fpout,"boundry location = %f \n",tmplayerptr->bl);
	  layernum++;
	  fprintf(fpout,"-----------------------------------------------\n\n");
	  tmplayerptr = tmplayerptr->nextptr;
	}

      fprintf(fpout,"-----------------------------------------------\n\n");
      fprintf(fpout,"Semi-infinite Layer: \n");
      fprintf(fpout, "\n Outer boundary is boundry %d \n",
	      (layernum - 1));
      
      ac = zsqrt(tmplayerptr->qn);
      fprintf(fpout,"pe = (%f %f), pm = (%f %f) \n",
	      tmplayerptr->per, tmplayerptr->pei,
	      tmplayerptr->pmr, tmplayerptr->pmi);
      fprintf(fpout,"layer number = %d \n", layernum);
      fprintf(fpout,"qn = (%f %f), yn = (%f %f) \n",
	      tmplayerptr->qn, tmplayerptr->yn);
      fprintf(fpout,"rn = (%f %f) \n\n", ac);
     fprintf(fpout,"-----------------------------------------------\n\n");

      /* describe first boundry condition */
      fprintf(fpout,"Boundary 1, kbc1 = %d, kbd1 = %d. \n",
	      modconptr->kbc1, modconptr->kbd1);
      fprintf(fpout,"Open boundry, semi-infinite adjacent layer \n");
      fprintf(fpout,"Principal branch specification for wx. Angle \n");
      fprintf(fpout,"apb = %f pi, vector vpb = (%f, %f). (wx*dot*vpb >= 20)\n",
	      modconptr->apb1, modconptr->vpb1);

      switch(modconptr->kbd1)
	{
	case 1:
	  fprintf(fpout,"Eigen condition: Inward wave solution only \n");
	  break;

	case 2:
	  fprintf(fpout,"Eigen condition: Outward wave solution only \n");
	  break;

	default:
	  printf("mode eigen condition not recognized \n");
	  exit(1);
	}
      fprintf(fpout,"Eigen condition- closed boundary. Outer layers ignored.");
      fprintf(fpout,"\n");
      fprintf(fpout,"Fixed surface admittance(TE)/impedance(TM)");
      fprintf(fpout," yb = (%f, %f) \n", modconptr->yb1);
      fprintf(fpout,"Fixed surface admittance(TE)/impedance(TM)");
      fprintf(fpout," zb = (%f, %f) \n \n", modconptr->zb1);

      /* describe second boundry condition */
      fprintf(fpout,"boundary %d, kbc2 = %d, kbd2 = %d. \n",
	      (layernum - 1), modconptr->kbc2, modconptr->kbd2);
      fprintf(fpout,"Open boundry, semi-infinite adjacent layer. \n");
      fprintf(fpout,"Principal branch specification for wx. Angle \n");
      fprintf(fpout,"apb = %f pi, vector vpb = (%f %f). (wx.vpb >= 20) \n",
	      modconptr->apb2, modconptr->vpb2);

      switch(modconptr->kbd2)
	{
	case 1:
	  fprintf(fpout,"Eigen condition: Inward wave solution only \n");
	  break;

	case 2:
	  fprintf(fpout,"Eigen condition: Outward wave solution only \n");
	  break;

	default:
	  printf("mode eigen condition not recognized \n");
	  exit(1);
	}
      fprintf(fpout," \n");
      fprintf(fpout,"Eigen condition- closed boundary. Outer layers ignored.");
      fprintf(fpout,"\n");
      fprintf(fpout,"Fixed surface admittance(TE)/impedance(TM)");
      fprintf(fpout," yb = (%f, %f) \n", modconptr->yb2);
      fprintf(fpout,"Fixed surface admittance(TE)/impedance(TM)");
      fprintf(fpout," zb = (%f, %f) \n\n", modconptr->zb2);

      fprintf(fpout,"end structure description \n");
      fprintf(fpout,"================================================== \n");
      fclose(fpout);
    }
  
  /* determine if it is a 1-D or a 2-D problem */
  if (caseptr->lxyopt == 0)
    {
      caseptr->kdof = 2;
      caseptr->kouf = MAX(caseptr->kouf, 3);
    }
  else
    {
      caseptr->kdofy = 2;
      caseptr->koufy = MAX(caseptr->koufy,3);
    }

  caseptr->ntheta = 1000;
  if (caseptr->dtheta != 0)
    {
      caseptr->ntheta = (int) ((caseptr->thetax/caseptr->dtheta) + 0.001);
    }

  if (caseptr->ntheta > 180)
    {
      printf("Too many points for far field %d \n",caseptr->ntheta);
      printf("Maximum 360 points permitted \n");
      exit(1);
    }
  return 0;
}
