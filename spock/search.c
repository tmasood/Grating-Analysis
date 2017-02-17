#include <stdio.h>
#include <math.h>
#include "util.h"
#include "case.h"
#include "layer.h"
#include "struct.h"
#include "modcon.h"
#include "output.h"

/* The set of guesses, at which the search for each root is started, */
/* is the single most important factor in determining which roots */
/* are found, and the order in which they are calculated.         */
/* The phase integral provides a basis for generating good guesses */
/* for complex modes of interest.                                  */

extern int charmatrix(struct LAYER *, struct STRUCT *, struct MODCON *, int);
extern int czerom(struct CASE *, struct STRUCT *, struct LAYER *,
		  struct OUTPUT *, struct MODCON *);

/* Search function calculates discrete modes of the waveguide. */
/* It calculates initial guesses for complex root search,  */
/* initializes necessary parameters and calls czerom for root */
/* search. For all the modes that are found it calculates */
/* propagation constants, admittances/impedances and phase */
/* integrals. It then prints summary of each mode parameters. */

int search(FILE *fpdb, struct CASE *caseptr, struct STRUCT *structptr,
           struct LAYER *headlayerptr, struct MODCON *modconptr,
	   struct OUTPUT *outputptr)
{
  struct LAYER *tmplayerptr;

  /* points representing guesses */
  double qr0;
  double qr1;
  double qr2;
  double qr3;
  double qr;

  double pr; /* total radian thickness */
  double aa;

  /* normalized phase integral */
  double pm0;
  double pm1;
  double pm2;
  double pm3;
  double pm;
  double qt;
  double d1, d2, a1, a2;

  dcomplex cc, ac, bc, dc;
  dcomplex CONE;

  int m, layernum;

  __real__ CONE = 1.0;
  __imag__ CONE = 0.0;

  /* Initial default values of root guesses. */
  /* Values of z in the neighborhood of root. */
  qr0 = 0.0;
  qr1 = caseptr->qznr;
  qr2 = 0.5;
  qr3 = 1.0; /* propagation constant of air */
  pr = 0.0;

  /* Minimum three layers required. Two outer boundries */
  /* and atleast one inner finite layer. */
  tmplayerptr = headlayerptr->nextptr;

  /* Keep qznr 9/10th of minimum qn and qr. */
  /* qr3 is the maximum value of qn (rn**2) over the inner layers */
  while (tmplayerptr->nextptr != NULL)
    {
      aa = __real__ tmplayerptr->qn;
      qr3 = MAX(aa, qr3);
      aa = 0.9 * aa;
      if ((aa >= 0.0) && (aa <= caseptr->qznr))
	{
	  caseptr->qznr = aa;
	}
      pr = pr + tmplayerptr->tr;
      tmplayerptr = tmplayerptr->nextptr;
    }

  /* q = 0 over the inner layers */
  __real__ cc = qr0;
  __imag__ cc = 0.0;
  structptr->qz = cc;
  charmatrix(headlayerptr, structptr, modconptr, 2);
  pm0 = __real__ structptr->phr / M_PI;

  /* qr1 has the lowest value of qn and qr over the inner layers */
  qr1 = caseptr->qznr;
  __real__ cc = qr1;
  __imag__ cc = 0.0;
  structptr->qz = cc;
  charmatrix(headlayerptr, structptr, modconptr, 2);
  pm1 = __real__ structptr->phr / M_PI;

  /* qr2 is mean value of qr1 and qr3 */
  qr2 = (qr1 + qr3) / 2.0;
  aa = __real__ (structptr->phr * structptr->phr);
  cc = pr * structptr->kc;
  qt = (__real__ cc * __real__ cc) + (__imag__ cc * __imag__ cc);

  if ((aa >= 0.001) && (qt >= 0.001))
    {
      qr2 = qr1 + aa / qt;
    }

  qr2 = (qr1 + qr2 + qr3) / 3.0;
  __real__ cc = qr2;
  __imag__ cc = 0.0;
  structptr->qz = cc;
  charmatrix(headlayerptr, structptr, modconptr, 2);
  pm2 = __real__ structptr->phr / M_PI;

  __real__ cc = qr3;
  __imag__ cc = 0.0;
  structptr->qz = cc;
  charmatrix(headlayerptr, structptr, modconptr, 2);
  pm3 = __real__ structptr->phr / M_PI;

  /* Prepare guesses and parameters for call to czerom. */
  switch (caseptr->kgss)
    {
    case 1:
      /* Use qzmr input value, or root values from previous */
      /* case as root guess. */
      for (m = 0; m < structptr->mn; m++)
	{
	  structptr->qz = caseptr->qzm[m];
	  charmatrix(headlayerptr, structptr, modconptr, 2);
	  caseptr->z0[m] = caseptr->qzm[m];
	  __real__ caseptr->phm[m] = __real__ structptr->phr / M_PI;
	  __imag__ caseptr->phm[m] = __imag__ structptr->phi / ED;
	  caseptr->fm[m] = -1.0;
	  caseptr->it[m] = 0;
	  caseptr->kr[m] = caseptr->km[m];
	}
      break;

    case 2:
      /* All initial guesses are the same; equal to max qr of inner finite */
      /* layers.                   */
      __real__ ac = qr3;
      __imag__ ac = caseptr->qzni;
      structptr->qz = ac;
      charmatrix(headlayerptr, structptr, modconptr, 2);
      __real__ bc = __real__ structptr->phr / M_PI;
      __imag__ bc = __imag__ structptr->phi / ED;

      for (m = 0; m < structptr->mn; m++)
	{
	  caseptr->qzm[m] = ac;
	  caseptr->phm[m] = bc;
	  caseptr->km[m] = caseptr->kgcz;
	  caseptr->z0[m] = ac;
	  caseptr->fm[m] = -1.0;
	  caseptr->it[m] = 0;
	  caseptr->kr[m] = caseptr->kgcz;
	}
      break;

    case 3:
      /* Individual guesses for qz. Simple quadratic approximation for */
      /* caseptr->qzm .vs. caseptr->phm for intended mode indices mm. */
      d1 = (qr1 - qr3)/ pm1;
      d2 = (qr2 - qr3)/ pm2;
      a2 = (d1 - d2) / (pm1 - pm2);
      if (a2 > 0.0)
	{
	  a2 = 0.0;
	}
      a1 = d1 - (a2 * pm1);
      /* Revise quadratic coefficients on basis of qr guess for avg mm. */
      pm = (((float) (structptr->mn))/ 2.0) + caseptr->pmfr;
      qr2 = qr3 + (a1 + (a2 * pm)) * pm;
      __real__ bc = qr2;
      __imag__ bc = caseptr->qzni;
      structptr->qz = bc;
      charmatrix(headlayerptr, structptr, modconptr, 2);
      pm2 = __real__ structptr->phr / M_PI;

      /* Revise coefficients for approximating quadratic */
      d2 = (qr2 - qr3) / pm2;
      a2 = (d1 - d2) / (pm1 - pm2);
      if (a2 > 0.0)
	{
	  a2 = 0.0;
	}
      a1 = d1 - (a2 * pm1);
      /* form initial guesses */
      for (m = 0; m < structptr->mn; m++)
	{
	  /* Expected phase integral pm = ms + po */
	  /* 's' is the mode spacing, often unity and po is the offset */
	  /* or bias. For simple guiding structures po = 1/4 to 1/2; but */
	  /* for a structure of three guides po = 3/2 may be more reasonable, */
	  /* and for the three lowest order modes 's' s may easily be much */
	  /* less than unity. The values of pm are used in the quadratic */
	  /* interpolation expression to generate the guesses */
	  pm = (caseptr->pmdm * ((float) m)) + caseptr->pmfr;
	  qr = qr3 + (a1 + (a2 * pm)) * pm;
	  __real__ ac = qr;
	  __imag__ ac = caseptr->qzni;
	  caseptr->qzm[m] = ac;
	  structptr->qz = ac;
	  charmatrix(headlayerptr, structptr, modconptr, 2);
	  __real__ caseptr->phm[m] = __real__ structptr->phr / M_PI;
	  __imag__ caseptr->phm[m] = __imag__ structptr->phi / ED;
	  caseptr->km[m] = caseptr->kgcz;
	  caseptr->z0[m] = ac;
	  caseptr->fm[m] = -1.0;
	  caseptr->it[m] = 0;
	  caseptr->kr[m] = caseptr->kgcz;
	}
      break;
      
    default:
      printf("kgss = %d not valid \n",caseptr->kgss);
      exit(1);
      break;
    }

  /* find root */
  czerom(caseptr, structptr, headlayerptr, outputptr, modconptr);

  for (m = 0; m < structptr->mn; m++)
    {
      caseptr->qzm[m] = caseptr->z0[m]; /* save z0 */
      caseptr->km[m] = caseptr->kr[m]; /* save kr into km */
      structptr->qz = caseptr->z0[m];
      charmatrix(headlayerptr, structptr, modconptr, 5);

      structptr->wzm[m] = structptr->wz; /* save wz */
      __real__ caseptr->phm[m] = __real__ structptr->phr/M_PI;
      __imag__ caseptr->phm[m] = __imag__ structptr->phi/ED;

      fprintf(fpdb," %f \t", caseptr->qzmr);
      fprintf(fpdb," %f \t", caseptr->qzmi);
      fprintf(fpdb," %f \t",__real__ caseptr->phm[m]);
      fprintf(fpdb,"(%f %g) \t",structptr->wz);
      if (caseptr->it[m] > caseptr->il)
	{
	  printf("Iteration limit exceeded \n");
	}
      if (caseptr->km[m] > 5)
	{
	  printf("\n \n \t\t **Root Search Converged** \n");
	  printf("   \t\t   --------------------- \n\n\n");
	  printf(" \t  phr \t \t \t phi \n");
	  printf(" \t  --- \t \t \t --- \n");
	  printf("(%f, %f) \t\t (%f, %f) \n\n",structptr->phr, structptr->phi);
	  printf("Renormalized Tranmission Matrix: \n ");
	  printf(" \t  sm0 \t \t \t sm1 \n");
	  printf(" \t  --- \t \t \t --- \n");
	  printf("(%f, %f) \t (%f, %f) \n\n",Smatrix.m0, Smatrix.m1);
	  printf(" \t  sm2 \t \t \t sm3 \n");
	  printf(" \t  --- \t \t \t --- \n");
	  printf("(%f, %f) \t (%f, %f)\n\n", Smatrix.m2, Smatrix.m3);
	  printf(" \t  dets \n");
	  printf(" \t  ---- \n");
	  printf("(%f, %f) \n\n",structptr->dets);
	  printf("Transmission Matrix: \n ");
	  printf(" \t  cm0 \t \t \t cm1 \n");
	  printf(" \t  --- \t \t \t --- \n");
	  printf("(%f, %f) \t (%f, %f) \n\n",Tmatrix.m0, Tmatrix.m1);
	  printf(" \t  cm2 \t \t \t cm3 \n");
	  printf(" \t  --- \t \t \t --- \n");
	  printf("(%f, %f) \t (%f, %f) \n\n",Tmatrix.m2, Tmatrix.m3);
	  dc = structptr->detc[0] - CONE;
	  printf("det = (%f, %f) \t detnorm = (%f, %f) \n\n",
		 dc, structptr->detc[1]);
	  /* for first semi-inifinite layer */
	  if (headlayerptr != NULL)
	    {
	      printf("First semi-infinite layer: \n");
	      printf("layer number = 1,  wx = (%f, %f),  yx = (%f, %f) \n\n",
		     headlayerptr->wx, headlayerptr->yx);
	    }

	  tmplayerptr = headlayerptr->nextptr;
	  /* wave parameters for each inner finite layers */
	  layernum = 2;
	  printf("Finite layers: \n");
	  while (tmplayerptr->nextptr != NULL)
	    {
	      printf("layer number = %d,  wx = (%f, %f),  yx = (%f, %f) \n",
		     layernum, tmplayerptr->wx, tmplayerptr->yx);
	      printf("ph = (%f, %f) \n", tmplayerptr->ph);
	      printf("cl = (%f, %f) (%f, %f) (%f, %f) \n",
		     tmplayerptr->cl0, tmplayerptr->cl1, tmplayerptr->cl2);
	      printf("cl3 = cl0 \n");
	      layernum++;
	      tmplayerptr = tmplayerptr->nextptr;
	    }
	  printf("\n");
	  /* for the last semi-infinite layer */
	  printf("Last semi-infinite layer: \n");
	  printf("layer number = %d, wx = (%f, %f), yx = (%f, %f) \n\n",
		 layernum, tmplayerptr->wx, tmplayerptr->yx);
	}
      else
	{
	  printf("Root search did not converge. \n");
	}
    }
  
  for (m = 0; m < structptr->mn; m++)
    {
      aa = __real__ caseptr->phm[m];
      if (0.0 > (aa + 0.5))
	{
	  caseptr->km[m] = caseptr->km[m] - 1;
	}
      if ((aa - 1.5) > (float)(structptr->mn - 1))
	{
	  caseptr->km[m] = caseptr->km[m] - 1;
	}

      fprintf(fpdb," \t %d \t",caseptr->km[m]);
      fprintf(fpdb," %d \n",caseptr->it[m]);
    }
  structptr->totallayers = layernum;
  return 0;
}
