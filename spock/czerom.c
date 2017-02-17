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
#include "output.h"
#include "modcon.h"

extern dcomplex eigeqf(struct STRUCT *, struct LAYER *,
		       struct CASE *, struct MODCON *, dcomplex, int);

/* finds complex roots of complex function eigeqf() using */
/* Muller-Traub method. */
/* eigeqf is the function whose roots are calculated */

int czerom(struct CASE *caseptr, struct STRUCT *structptr,
           struct LAYER *headlayerptr, struct OUTPUT *outputptr,
	   struct MODCON *modconptr)
{
  int n, nl, ii, nt, nr, ks;
  int l, kant, kout;
  int i, j, k;
  float epsz, epsf, epqz, epqf, bb, aa1;
  double rr, ff, fqz;
  double fq[4], *fqr;
  
  dcomplex *zz, cc, cg, rotc, cz, zg, *fz, tmpconst;
  dcomplex zi[4], cd, fi[4], fcz, dz;
  dcomplex ca, cb, aa, dfab, dfbc, ddfc, dfdz;
  dcomplex *fa, *fb, *fc;
  dcomplex *za, *zb, *zc;

  dcomplex CONE, CZERO;
  
  __real__ CONE = 1.0;
  __imag__ CONE = 0.0;

  __real__ CZERO = 0.0;
  __imag__ CZERO = 0.0;

  /* Complex rotation factor R with an angle of about 120 degrees. */
  /* R = -0.6 + i0.8 */
  __real__ rotc = -0.6;
  __imag__ rotc = 0.8;

  fa = &fi[0];
  fb = &fi[1];
  fc = &fi[2];
  fz = &fi[3];
  za = &zi[0];
  zb = &zi[1];
  zc = &zi[2];
  zz = &zi[3];
  fqr = &fq[3];

  /* nl is the number of roots to be searched */
  nl = structptr->mn;
  if (nl == 0)
    {
      printf("No root searched \n");
      exit(1);
    }

  /* set lower limit on epsilon */
  aa1 = 100.0 * AM;
  epsz = caseptr->eps1;
  epsf = caseptr->eps2;
  if (epsz < aa1)
    {
      epsz = aa1;
    }
  if (epsf < aa1)
    {
      epsf = aa1;
    }
  ii = caseptr->il;

  /* squared epsilon form used in convergence tests */
  epqz = epsz * epsz;
  epqf = epsf * epsf;
  cg = rotc;
  kout = caseptr->ktuo;
  kant = 0;
  
  if (caseptr->lxyopt <= 0)
    {
      if ((kout >= 1) && (outputptr->modout > 0))
	{
	  printf("Subroutine CZEROM iterations \n");
	}
      if ((kout >= 3) && (outputptr->modout > 0))
	{
	  printf("n \t kr \t i \ndz[i] \t z[i] \nf[z] \t");
	  printf("f[z]-reduced \t f-sqd \n\n");
	}
    }

  /* main loop searching for nl roots */
  for (n = 0; n < nl; n++)
    {
      nt = nr = n;
      ks = caseptr->kr[n];
      *zz = caseptr->z0[n];

      if ((ks < 0) || (ks > 7))
	{
	  /* no search for fixed root. Just set parameters */
	  caseptr->it[n] = 0;
	  dz = CZERO;
	  fcz = eigeqf(structptr, headlayerptr, caseptr, modconptr, *zz, nr);
	  *fz = -1.0 * CZERO;
	  *fqr = -0.0;
	  /* bypass root search */
	}
      else 
	{
	  /* proceed with root search */
	  /* setup guesses. Three on a circle about z0, and one near z0 */
	  i = 0;
	  if (ks == 0)
	    {
	      *zz = CONE;
	    }
	  /* cc is the complex conjugate of z0 */
	  cc = zcnjg(*zz);
	  /* set a max radius for circle of guesses about an origin */
	  /* ro = |z0| */
	  rr = zabs(*zz);
	  /* rmax = (1 / (1 + 2r0)) */
	  rr = 0.5 / (rr + 0.5);

	  /* It is important to be able to control the distance of the */
	  /* three points about the initial guess. It is convenient to */
	  /* choose a radius which is negative power of 10, say 10**(-k). */
	  /* Where k is chosen to reflect the quality of the guess. A poor */
	  /* guess may use a k of about 0 to 2, while an excellent guess */
	  /* may justify a k of 3 to 5. The radius should be larger than the */
	  /* convergence criteria being used, and it should be smaller */
	  /* than the expected spacing between roots. */
	  /* Use a fraction of max radius, powers of (0.1)**caseptr->kr[n] */
	  ff = 1.0;
	  if (ks > 1)
	    {
	      ff = pow(0.1,ks);
	    }

	  /* Set radius vector for first guess. */
	  /* wj = R * rw,  rw = 10**k * rmax, j = 1,2,3 */
	  cz = ff * rr * cg;

	  /* ****************************************************** */
	  /* The three points wj are transformed to become points about */
	  /* z0 by the Mobius transformation */
	  /* Generate the four guesses */
	  for (k = 0; k < 4; k++)
	    {
	      /* Unitary Mobius transformation to points about z0 */
	      zg = (*zz + cz)/(CONE - (cc * cz));
	      zi[k] = zg;
	      fcz = eigeqf(structptr, headlayerptr, caseptr, modconptr,
			   zg, nr);

	      /* factor out known and fixed roots, but not other */
	      /* roots that are poorly known or current root. */
	      /* Form denominator polynomial */
	      cd = CONE;
	      for (l = 0; l < nl; l++)
		{
		  if ((caseptr->kr[l] <= 3) || (l == n))
		    {
		      break;
		    }
		  else
		    {
		      ca = zg - caseptr->z0[l];
		      /* if too near a known root accept as a */
		      /* multiple root. Avoids inflating reduced */
		      /* function or numerical indeterminacy. */
		      /* Chordal metric on Riemann Sphere for z-plane */
		      bb = (1.0 + zabsq(zg)) * (1.0 + zabsq(caseptr->z0[l]));
		      if (zabsq(ca) > (epqz * bb))
			{
			  cd = cd * ca;
			}
		      else
			{
			  *zz = zg;
			  *fz = fcz;
			  *fqr = fq[k];
			  ks = caseptr->kr[l];
			  caseptr->it[n] = 0;
			  goto L500;
			}
		    }
		}

	      /* divide out polynomial of known roots */
	      cd = fcz / cd; 
	      fi[k] = cd;
	      fq[k] = zabsq(cd);

	      /* Rotate radius vector for next guess. */
	      /* For last guess use small radius for */
	      /* approx zero avg of the four guesses */
	      cz = cz * rotc;
	      if (k >= 2)
		{
		  __real__ tmpconst = -0.019;
		  __imag__ tmpconst = 0.064;
		  cz = tmpconst * cz;
		}
	      if (caseptr->lxyopt <= 0)
		{
		  /* output each initial guess iterate */
		  if ((kout >= 4) && (outputptr->modout > 0))
		    {
		      printf("n \t ks \t i \t C0 \t zi \n");
		      printf("%d \t %d \t %d \t (%f %f) \t (%f %f) \n",
			     n, ks, i, CZERO, zi[k]);
		      printf("fcz \t fi \t fq \n");
		      printf("(%f %f) \t (%f %f) \t %f \n",
			     fcz, fi[k], fq[k]);
		    }
		}
	    } /* for (k=0; k<4; k++):- transformation of points */

	  /* ****************************************************** */
	  /* set radius vector for next guesses */
	  cg = -cg * rotc;
	  /* guesses completed */
	  /* Only if the three iterates are in the neighborhood of an */
	  /* isolated root will the successive iterates always be smaller */
	  /* and the convergence order of 1.84 be realized.           */
	  /*                                                          */
	  /* reorder guesses in decreasing absq of reduced function */
	  /* fq0 >= fq1 >= fq2 >= fq3 */
	  for (i = 0; i < 3; i++)
	    {
	      j = i + 1;
	      for (k = j; k < 4; k++)
		{
		  if (fq[k] >= fq[i])
		    {
		      ca = zi[i];
		      cb = fi[i];
		      aa = fq[i];
		      zi[i] = zi[k];
		      fi[i] = fi[k];
		      fq[i] = fq[k];
		      zi[k] = ca;
		      fi[k] = cb;
		      fq[k] = aa;
		    }
		}
	    }
	  /* The largest iterate is discarded */
	  /* dump guess with largest fq and push the rest up */
	  for (i = 0; i < 3; i++)
	    {
	      j = i + 1;
	      zi[i] = zi[j];
	      fi[i] = fi[j];
	      fq[i] = fq[j];
	    }

	  /* guesses ordered, used as initial iterates */
	  /* main muller iteration loop */
	  ks = 4;
	  caseptr->kr[n] = ks;
	  for (i = 0; i < ii; i++)
	    {
	      /* calculate next muller iterate. Traub formulas */
	      /* Iteration based on reduced function. */
	      /* fa is the largest and fc the smallest */
	      dfab = (*fa - *fb) / (*za - *zb);
	      dfbc = (*fb - *fc) / (*zb - *zc);
	      ddfc = (dfab - dfbc) / (*za - *zc);
	      cb = dfbc + (*zc - *zb) * ddfc;
	      ca = zsqrt(((cb * cb) - (4.0 * *fc * ddfc)));
	      cc = cb + ca;
	      cd = cb - ca;
	      /* select largest denominator */
	      if (zabsq(cc) < zabsq(cd))
		{
		  cc = cd;
		}
	      dfdz = 0.5 * cc;
	      dz = -(*fc)/dfdz;
	      /* new iterate approximation for root */
	      *zz = *zc + dz;
	      /* evaluate function at new iterate */
	      fcz = eigeqf(structptr, headlayerptr, caseptr, modconptr,
			   *zz, nr);
	      *fz = fcz;
	      fqz = zabsq(fcz);

	      /* factor out known roots from polynomial denominator */
	      cd = CONE;
	      for (l = 0; l < nl; l++)
		{
	      /* but not guesses that are poorly known, or is */
	      /* current root */
		  if ((caseptr->kr[l] <= 3) || (l == n))
		    {
		      break;
		    }
		  else
		    {
		      ca = *zz - caseptr->z0[l];
		      /* if too near a known root accept as multiple root */
		      /* avoids inflating reduced function or numerical */
		      /* indeterminacy. Chordal metric on Riemann sphere */
		      /* for z-plane */
		      if (caseptr->kr[l] > 4)
			{
			  bb = bb = (1.0 + zabsq(*zz)) * 
			    (1.0 + zabsq(caseptr->z0[l]));
			  if (zabsq(ca) <= (epqz * bb))
			    {
			      ks = caseptr->kr[l];
			      caseptr->it[n] = i;
			      goto L500;
			    }
			  else
			    {
			      cd = cd * ca;
			    }
			}
		      else
			{
			  cd = cd * ca;
			}
		    }
		}  /* for (l = 0; l < nl; l++) */
	      
	      /* divide out polynomials of known roots */
	      *fz = *fz / cd;
	      *fqr = zabsq(*fz);
	      printf("fqr = %g fqz = %g i = %d epqf = %g \n",*fqr, fqz, i, epqf);
	      /* test for convergence */
	      /* convergence test on magnitude of reduced function */
	      if (*fqr < epqf)
		{
		  ks = 5;
		}
	      /* convergence test on magnitude of actual function */
	      if (fqz < epqf) 
		{
		  ks = 6;
		}
	      /* convergence test on iterate increment dz */
	      /* chordal metric on Riemann sphere for z-plane */
	      if (zabsq(dz) <= (epqz*(1.0 + (zabsq(*zz)*zabsq(*zz)))))
		{
		  ks = 7;
		}

	      if (ks < 5)
		{
		  /* not converged. Complete iteration. Update array */
		  /* of iterates. Dump largest old iterate. Insert new */
		  /* iterate in order */
		  for (j = 0; j < 2; j++)
		    {
		      k = j + 1;
		      if (*fqr <= fq[k])
			{
			  zi[j] = zi[k];
			  fi[j] = fi[k];
			  fq[j] = fq[k];
			}
		      else
			{
			  break;
			}
		    }
		  zi[j] = *zz;
		  fi[j] = *fz;
		  fq[j] = *fqr;
		  if (caseptr->lxyopt <= 0)
		    {
		      /* output each iteration */
		      if ((kout >= 4) && (outputptr->modout > 0))
			{
			  printf("%d \t %d \t %d \t (%f %f) \t (%f %f) \n",
				 n, ks, i, dz, zz);
		      printf("(%f %f) \t (%f %f) \t %f \n",
			     fcz, *fz, *fqr);
			}
		      /* each iteration completed */
		    }
		}
	      else
		{
		  break;
		}
	    } /* for (i = 0; i < ii; i++) */

	  /* iteration limit */
	  caseptr->it[n] = i;
	  kant = kant + 1;

L500:
	  caseptr->z0[n] = *zz;
	  caseptr->fm[n] = sqrt(fqz);
	  caseptr->kr[n] = ks;
	  if (ks >= 5)
	    {
	      kant = 0;
	    }
	} /* end else block:- if ((ks < 0) || (ks > 7) */

      /* ****************************************************** */
      if (caseptr->lxyopt > 0)
	{
	  if (kant != 0)
	    {
	      printf("Root failed to converge for mode. \n");
	      printf("%d \t %d \t %d \n",n, caseptr->it[n], ks);
	      if (n == 0)
		{
		  exit(1);
		}
	    }
	}
      else
	{
	  if ((kout >= 3) && (outputptr->modout > 0))
	    {
	      printf("%d \t %d \t %d \n",n, caseptr->kr[n], caseptr->it[n]);
	      printf("(%g, %g) \t (%g, %g) \n",dz, caseptr->z0[n]);
	      printf("(%g, %g) \t (%g, %g) \t %g \n\n",fcz, *fz, *fqr);
	      printf("***************************************\n\n");
	    }
	  if (kant >= 2)
	    {
	      /* break out of top for loop:- for(n=0;n<nl;n++) */
	      /* The main loop for root search */
	      break;
	    }
	}
    } /* main search loop for nl roots:- for(n=0;n<nl;n++) */

  /* ************************************************************* */
  return 0;
}
