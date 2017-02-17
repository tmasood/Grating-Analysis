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
#include "local_complex.h"

extern int laymat(double, dcomplex, dcomplex, dcomplex,
		  dcomplex *, dcomplex *, dcomplex *);

extern int ftfeld(char *, double *, dcomplex *, dcomplex *, int,
                  struct STRUCT *, struct OUTPUT *, struct CASE *,
		  struct LAYER *);

int fields(struct CASE *caseptr, struct STRUCT *structptr,
           struct LAYER *headlayerptr, struct MODCON *modconptr,
           struct OUTPUT *outputptr, char *filename)
{
  dcomplex kc0;
  struct LAYER *tmplayerptr, *tmp1layerptr,  *lastlayerptr, *tmpblayerptr;
  dcomplex fy1, fy2, fz1, fz2;
  dcomplex fg1, fg2, gf1, gf2;
  dcomplex dc, ac, bc, cc;
  dcomplex fyz, fxx, fyx, fzx;
  dcomplex pxx, pzx, tm1, tm2, tm3;
  dcomplex fyxft[MAXPOINTS], pzxft[MAXPOINTS];

  double fn, gn, aa, dx, tx, xx;
  double dp, x2, delx, imax, iden;
  double delw, dl, x1, infx, xxft[MAXPOINTS];
  int tmp;
  int kbfg, m, icnt, i, ix, jl;
  int kdo0;

  kc0 = structptr->kc * structptr->k0;
  iden = 1000.0;
  imax = 15;
  infx = 3.0;
  delw = 1.0e-4;
  delx = 1.01e-4;

  /* Calculate tangential fields and transverse poynting power. */
  for (m = 0; m < structptr->mn; m++)
    {
      structptr->mm = structptr->mo + m;
      structptr->qz = caseptr->qzm[m];

      /* calculate propagation and transformation matrix variables */
      kdo0 = 5;
      charmatrix(headlayerptr, structptr, modconptr, kdo0);
      __real__ caseptr->phm[m] = __real__ structptr->phr / M_PI;
      __imag__ caseptr->phm[m] = __imag__ structptr->phi / ED;

      /* Skip field calculations for nonconverged eigen or poor qz */
      /* value. */
      if (caseptr->km[m] < caseptr->kmdo)
	{
	  printf("Field of nonconverged neff is not calculated \n");
	  exit(1);
	}
      if (caseptr->kouf >= 3)
	{
	  if (caseptr->lxyopt <= 0)
	    {
	      aa = 0.5 * caseptr->aqz;
	    }
	}
      /* setup two provisional sets of tangential boundary fields. */
      /* from outer boundary fields and layer matrices saved in fy, */
      /* fz from first (top) boundary transform forward to second */
      /* (bottom) boundary .*/
      /* ... and matrices saved in gy, gz from second (bottom) boundary */
      /* transform backward to first (top) boundary */

      /* set outer boundary fields */
      tmplayerptr = headlayerptr;
      fy1 = tmplayerptr->fy;
      fz1 = tmplayerptr->fz;

      while (tmplayerptr->nextptr != NULL)
	{
	  tmplayerptr = tmplayerptr->nextptr; 
	}
      lastlayerptr = tmplayerptr;
      /* tangential fields at the second boundary */
      fy2 = tmplayerptr->fy;
      fz2 = tmplayerptr->fz;
      tmplayerptr->gy = fy2;
      tmplayerptr->gz = fz2;
      /* accumulate normalization */
      fn = zabsq(fy1) + zabsq(fz1);
      gn = zabsq(fy2) + zabsq(fz2);

      /* no finite layer transformation for single boundry case */
      /* calculate matrix transformation across finite layers in */
      /* both directions */
      if (headlayerptr->nextptr != NULL)
	{
	  /* Matrix transformation across finite layers, both */
	  /* directions */
	  tmplayerptr = headlayerptr;
	  tmpblayerptr = lastlayerptr;
	  while (tmplayerptr->nextptr->nextptr != NULL)
	    {
	      tmplayerptr->nextptr->fy = ((tmplayerptr->nextptr->cl0 * 
					   tmplayerptr->fy) +
					  (tmplayerptr->nextptr->cl2 *
					   tmplayerptr->fz));
	      tmplayerptr->nextptr->fz = ((tmplayerptr->nextptr->cl1 * 
					   tmplayerptr->fy) +
					  (tmplayerptr->nextptr->cl0 *
					   tmplayerptr->fz));
	      tmpblayerptr->prevptr->gy = ((tmpblayerptr->prevptr->cl0 * 
					    tmpblayerptr->gy) -
					   (tmpblayerptr->prevptr->cl2 *
					    tmpblayerptr->gz));
	      tmpblayerptr->prevptr->gz = (-(tmpblayerptr->prevptr->cl1 * 
					     tmpblayerptr->gy) +
					   (tmpblayerptr->prevptr->cl0 *
					    tmpblayerptr->gz));

	      /* accumulate normalization */
	      fn = fn + zabsq(tmplayerptr->nextptr->fy) + 
		zabsq(tmplayerptr->nextptr->fz);
	      gn = gn + zabsq(tmpblayerptr->prevptr->gy) + 
		zabsq(tmpblayerptr->prevptr->gz);
	      tmplayerptr = tmplayerptr->nextptr;
	      tmpblayerptr = tmpblayerptr->prevptr;
	    }
	}

      /* Fields are normalized to RMS (magnitude) of boundary values. */
      /* This is an arbitrary choice, other normalizations are */
      /* possible. */
      gn = sqrt(2.0 * (float)structptr->totallayers / gn);
      fn = sqrt(2.0 * (float)structptr->totallayers / fn);

      tmplayerptr = headlayerptr;
      tmpblayerptr = lastlayerptr;
      while(tmplayerptr->nextptr != NULL)
	{
	  tmplayerptr->fy = fn * tmplayerptr->fy;
	  tmplayerptr->fz = fn * tmplayerptr->fz;
	  tmpblayerptr->gy = gn * tmpblayerptr->gy;
	  tmpblayerptr->gz = gn * tmpblayerptr->gz;

	  tmplayerptr = tmplayerptr->nextptr;
	  tmpblayerptr = tmpblayerptr->prevptr;
	}
      kbfg = 1;
      
      /* f * g cross products for reciprocity and eigen checks */
      tmplayerptr = headlayerptr;
      tmpblayerptr = lastlayerptr;
      fg1 = tmplayerptr->fy * tmplayerptr->nextptr->gz;
      fg2 = tmpblayerptr->prevptr->fy * tmpblayerptr->gz;
      gf1 = tmplayerptr->nextptr->gy * tmplayerptr->fz;
      gf2 = tmpblayerptr->gy * tmpblayerptr->prevptr->fz;

      if (caseptr->lxyopt <= 0)
	{
	  if (caseptr->kouf >= 2)
	    {
	      /* Output summary of provisional field sets f */
	      /* and g. Reciprocity and eigen checks using */
	      /* poynting cross products f * g. */
	      if (caseptr->kouf >= 3)
		{
		  dc = fg2 - gf2;
		  dc = fg1 - gf1;
		  ac = fg2 - fg1;
		  bc = gf2 - gf1;
		  cc = ac - bc;
		  /* Reciprocal transmission and reflection coefficients.*/
		  /* Not yet implemented. */
		  if (caseptr->kouf >= 4)
		    {
		      /* Output both f and g sets of tangential boundary */
		      /* fields. */
		      tmplayerptr = headlayerptr;
		      while (tmplayerptr->nextptr != NULL)
			{
			  dc = (tmplayerptr->fy * tmplayerptr->nextptr->gz) -
			    (tmplayerptr->nextptr->gy * tmplayerptr->fz);
			  tmplayerptr = tmplayerptr->nextptr;
			}
		    }
		}
	    }
	}
      /* Use both, one or the other, or the average of the */
      /* two field solutions */
      /* Also calculate time-avg transverse poynting power */
      if (modconptr->kbc0 >= 3)
	{
	  /* A true eigen-mode case, f and g solutions equivalent */
	  /* Form average of f and g sets as eigen-function solutions */
	  /* Standardize gy solution to be the same as fy at first */
	  /* (top) layer */
	  tmplayerptr = headlayerptr;
	  cc = tmplayerptr->fy / tmplayerptr->nextptr->gy;
	  while (tmplayerptr->nextptr != NULL)
	    {
	      tmplayerptr->fy = 0.5 * (tmplayerptr->fy +
				       (cc * tmplayerptr->nextptr->gy));
	      tmplayerptr->fz = 0.5 * (tmplayerptr->fz +
				       (cc * tmplayerptr->nextptr->gz));
	      tmplayerptr->px = tmplayerptr->fy * zcnjg(tmplayerptr->fz);
	      tmplayerptr = tmplayerptr->nextptr;
	    }
	  kbfg = 3;

	}
      else if (modconptr->kbc0 == 2)
	{
	      /* Exchange f and g sets, for use with second field set */
	      /* gy, gz. A return entry point for using g set after */
	      /* using f set, kbc0 = 0 */
	      tmplayerptr = headlayerptr;
	      while (tmplayerptr->nextptr != NULL)
		{
		  ac = tmplayerptr->nextptr->gy;
		  bc = tmplayerptr->nextptr->gz;
		  tmplayerptr->nextptr->gy = tmplayerptr->fy;
		  tmplayerptr->nextptr->gz = tmplayerptr->fz;
		  tmplayerptr->fy = ac;
		  tmplayerptr->fz = bc;
		  tmplayerptr->px = ac * zcnjg(bc);
		  tmplayerptr = tmplayerptr->nextptr;
		}
	      kbfg = 2;
	}

      else
	{
	  /* Use first field. Set fy, fz. Only set if kbc0 = 1 */
	  tmplayerptr = headlayerptr;
	  while (tmplayerptr->nextptr != NULL)
	    {
	      tmplayerptr->px = tmplayerptr->fy * 
		zcnjg(tmplayerptr->fz);
	      tmplayerptr = tmplayerptr->nextptr;
	    }
	}

      /* save the field solution qz for the initial guess in the */
      /* next loop */
      caseptr->qzm[m] = structptr->qz;
      
      if (caseptr->lxyopt > 0)
	{
	  tmplayerptr = headlayerptr;
	  while (tmplayerptr->nextptr != NULL)
	    {
	      tmplayerptr->fym[m] = tmplayerptr->fy;
	      tmplayerptr->fzm[m] = tmplayerptr->fz;
	      tmplayerptr = tmplayerptr->nextptr;
	    }
	  return (0);
	}
      /* summary output description of type of field solution */
      /* and listing of tangential fields at boundaries. */

      if (kbfg == 2)
	{
	  tmplayerptr = headlayerptr;
	  tmpblayerptr = lastlayerptr->prevptr;
	  cc = (fz1 * tmplayerptr->fy) - (tmplayerptr->fz * fy1);
	  dc = (fz2 * tmpblayerptr->fy) - (tmpblayerptr->fz * fy2);
	}
	  
      /* output tangential boundary fields and transverse poynting */
      /* average power. */
      /* time average net poynting power into structure from the */
      /* outside. */
      /* output wave amplitudes at outer boundaries, transmission */
      /* and reflection coefficient not yet implemented */
      /* Calculate power and energy relations for field solution */
      powers();

      if ((caseptr->kdof <= 1) || (caseptr->kouf <= 1))
	{
	  tmplayerptr = headlayerptr;
	  while (tmplayerptr->nextptr != NULL)
	    {
	      tmplayerptr->fym[m] = tmplayerptr->fy;
	      tmplayerptr->fzm[m] = tmplayerptr->fz;
	      tmplayerptr = tmplayerptr->nextptr;
	    }
	  return (0);
	}

      /* calculate and output fields as a function of x. */
      /* Start loop over all layers. */
      /* Initialize counter for field fourier transformation array */
      icnt = 0;
      if (caseptr->dxin != 0.0)
	{
	  dx = caseptr->dxin;
	}

      tmplayerptr = headlayerptr;
      while (tmplayerptr != NULL)
	{
	  if (tmplayerptr == headlayerptr)
	    {
	      if (modconptr->kbc1 >= 2)
		{
		  tmp1layerptr = headlayerptr;
		  while (tmp1layerptr->nextptr != NULL)
		    {
		      tmp1layerptr->fym[m] = tmp1layerptr->fy;
		      tmp1layerptr->fzm[m] = tmp1layerptr->fz;
		      tmp1layerptr = tmp1layerptr->nextptr;
		    }
		  return (0);
		}
	      /* For first semi infinite layer if present */
	      /* set up an initial boundary and fields at */
	      /* infinity x. factors of e or radians into semi */
	      /* infinite layer */
	      tx = infx / (delw + zabs(tmplayerptr->wx * kc0));
	      x1 = tmplayerptr->xl - tx;
	      laymat(structptr->k0*tx, tmplayerptr->wx, structptr->kc,
		     tmplayerptr->yn, &tm1, &tm2, &tm3);
	      fyx = (tm1 * tmplayerptr->fy) - (tm3 * tmplayerptr->fz);
	      fzx = -(tm2 * tmplayerptr->fy) + (tm1 * tmplayerptr->fz);
	    }
	  /* for other layers store the fields */
	  else
	    {
	      x1 = tmplayerptr->prevptr->xl;
	      fyx = tmplayerptr->prevptr->fy;
	      fzx = tmplayerptr->prevptr->fz;
	    }

	  if (tmplayerptr == lastlayerptr)
	    {
	      if (modconptr->kbc2 >= 2)
		{
		  tmp1layerptr = headlayerptr;
		  while (tmp1layerptr->nextptr != NULL)
		    {
		      tmp1layerptr->fym[m] = tmp1layerptr->fy;
		      tmp1layerptr->fzm[m] = tmp1layerptr->fz;
		      tmp1layerptr = tmp1layerptr->nextptr;
		    }
		  return (0);
		}
	      /* for second semi-infinite layer if present. */
	      /* setup a final boundary at infx factors */
	      /* of e or radians into semi-infinite layer */
	      tx = infx / (delw + zabs(tmplayerptr->wx * kc0));
	      x2 = tmplayerptr->prevptr->xl + tx;
	    }
	  else
	    {
	      /* setup final x as xl, except the semi-infinite layer */
	      x2 = tmplayerptr->xl;
	    }

	  /* skip over logic to compute and adjust dx. dx is */
	  /* input as dxin */
	  if (caseptr->dxin == 0.0)
	    {
	      tx = x2 - x1;
	      dx = 1.0/(iden * (delw + zabs(tmplayerptr->wx * kc0)));
	      dx = MAX3(dx, (tx/(double) imax), delx);

	      /* Round-off dx to a two digit factor times a power */
	      /* of ten. Rounded values of 1.0, 1.5, 2.0, 2.5, 4.0, */
	      /* 5.0, 7.5, 10.0. About 1/3 rounds up and 2/3 rounds */
	      /* down. Log10 of dx. Integer and fractional part */
	      dl = log10(dx);
	      jl = (int) dl;
	      dl = dl - (double) jl;
	      /* check that fractional part is negative */
	      if (dl <= 0.0)
		{
		  dl = dl + 1.0;
		  jl = jl - 1;
		}

	      /* power of 10 factor, and round the fractional */
	      /* part */
	      dp = pow(10.0,(double) jl);
	      dx = 1.0;
	      if (dl >= 0.96)
		{
		  dx = 10.0;
		}
	      else if (dl >= 0.81)
		{
		  dx = 7.5;
		}
	      else if (dl >= 0.67)
		{
		  dx = 5.0;
		}
	      else if (dl >= 0.53)
		{
		  dx = 4.0;
		}
	      else if (dl >= 0.36)
		{
		  dx = 2.5;
		}
	      else if (dl >= 0.26)
		{
		  dx = 2.0;
		}
	      else if (dl >= 0.11)
		{
		  dx = 1.5;
		}
	      dx = dx * dp;
	    }
	  /* round down initial x1 to integer number of dx */
	  /* and set up fields */
	  tx = x1;
	  /* reset x1 to the nearest integral multiple of dx */
	  /* not exceeding x1 in magnitude */
	  x1 = dx * (int) ((x1+0.000001)/dx);
	  if ((tx+0.000001) < x1)
	    {
	      x1 = x1 - dx;
	    }
	  tx = x1 - tx;
	  if (tx != 0.0)
	    {
	      /* compute transformation matrix for thickness */
	      /* tx to backup to rounded value of x1. */
	      laymat(structptr->k0*tx, tmplayerptr->wx, structptr->kc,
		     tmplayerptr->yn, &tm1, &tm2, &tm3);
	      dc = tm1*fyx + tm3*fzx;
	      fzx = tm2*fyx + tm1*fzx;
	      fyx = dc;
	      pxx = fyx * zcnjg(fzx);
	      fxx = tmplayerptr->yz * fyz;
	      pzx = zcnjg(fxx) * fyx;
	    }

	  /* transformation matrix for dx thickness in layer L */
	  laymat(structptr->k0*dx, tmplayerptr->wx, structptr->kc,
		 tmplayerptr->yn, &tm1, &tm2, &tm3);

	  /* calculate and output fields from x1 to x2, increments */
	  /* of dx */
	  ix = (int)(((x2 - x1)/dx) + 0.5);
	  /* start do-loop incrementing dx through each layer */
	  for (i=1; i<=ix; i++)
	    {
	      xx = x1 + ((double) i * dx);
	      if (xx > (x2+0.000001))
		{
		  break;
		}
	      dc = tm1*fyx + tm3*fzx;
	      fzx = tm2*fyx + tm1*fzx;
	      fyx = dc;
	      pxx = fyx * zcnjg(fzx);
	      /* calculate and output fx and pz */
	      fxx = structptr->k0 * tmplayerptr->yz * fyx;
	      pzx = 0.5 * zcnjg(fxx) * fyx;
	      xxft[icnt] = xx;
	      fyxft[icnt] = fyx;
	      pzxft[icnt] = pzx;

	      /* now recompute transformation matrix using the */
	      /* new value of dx. */
	      /* dx = dx + dx; */
	      /* laymat(structptr->k0*dx, tmplayerptr->wx, structptr->kc,
		 tmplayerptr->yn, &tm1, &tm2, &tm3); */
	      icnt++;	
	    }

	  tmplayerptr = tmplayerptr->nextptr;
	}
      if ((caseptr->nfplt > 0) || (caseptr->ffplt > 0))
	{
	  ftfeld(filename, xxft, fyxft, pzxft, icnt, structptr, outputptr, caseptr, headlayerptr);
	}
    }
  return (0);
}
