#include <stdio.h>
#include <math.h>
#include "util.h"
#include "struct.h"
#include "layer.h"
#include "case.h"
#include "modcon.h"
#include "charmatrix.h"

extern int charmatrix(struct LAYER *, struct STRUCT *, struct MODCON *, int);

/* Calculate value of eigen equation function from system */
/* characteristic matrices and boundry conditions. */
dcomplex eigeqf(struct STRUCT *structptr, struct LAYER *headlayerptr,
		struct CASE *caseptr, struct MODCON *modconptr, 
		dcomplex zz, int nr)
{
  dcomplex qz, kc, eige, ac;
  double phmm, wf;
  int k;

  structptr->qz = zz;
  structptr->mm = structptr->mo + nr;
  charmatrix(headlayerptr, structptr, modconptr, 4);
  phmm = __real__ structptr->phr / M_PI;

  /* select eigen equation function from S-Matrix elements */
  eige = Smatrix.m3;
  wf = 1.0;

  if (caseptr->keif > 1)
    {
      /* modify eigeqf by a renormalization relative to unity. */
      /* divide by rms magnitudes of the 4 elements of S-matrix */
      wf = 0.0;
      for (k = 0; k < 4; k++)
	{
	  switch (k)
	    {
	    case 0:
	      ac = Smatrix.m0;
	      break;
	    case 1:
	      ac = Smatrix.m1;
	      break;
	    case 2:
	      ac = Smatrix.m2;
	      break;
	    case 3:
	      ac = Smatrix.m3;
	      break;
	    }
	  wf = wf + (__real__ ac * __real__ ac) + (__imag__ ac * __imag__ ac);
	}
      wf = sqrt(4.0/wf);
      eige = wf * eige;
      if (caseptr->keif > 2)
	{
	  /* modify eigeqf to bias search toward intended mode index. */
	  /* Multiply by positive real weight function. Square-root */
	  /* difference between phase integral and intended mode index */
	}
    }
  return (eige);
}
