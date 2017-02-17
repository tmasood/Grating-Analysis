#include <stdio.h>
#include <math.h>
#include "util.h"
#include "case.h"
#include "layer.h"
#include "struct.h"
#include "modcon.h"
#include "output.h"

/* calculate the full width half max intensity of the */
/* near and far fields.                               */

int width(int n, double *x, double *fy, double *fwhp)
{

  dcomplex t;
  double s[MAXPOINTS], q;
  double xhpl, xhpu;
  double fy2[MAXPOINTS];
  int i, ill, ilu;
  
  /* normalize the field */
  q = 0.0;
  for (i = 0; i < n; i++)
    {
      q = MAX(q,fy[i]);
    }

  for (i = 0; i < n; i++)
    {
      fy2[i] = fy[i]/q;
    }

  /* find the lower half power point */
  for (i = 0; i < n; i++)
    {
      if (fy2[i] > 0.5)
	{
	  break;
	}
    }
  ill = i - 1;
  ilu = i;
  xhpl = x[ill] + (x[ilu] - x[ill])*(0.5 - fy2[ill])/(fy2[ilu] - fy2[ill]);

  /* find the upper half power point */
  for (i = ilu; i < n; i++)
    {
      if (fy2[i] < 0.5)
	{
	  break;
	}
    }
  ill = i - 1;
  ilu = i;
  xhpu = x[ilu] + (x[ilu] - x[ill])*(0.5 - fy2[ilu])/(fy2[ilu] - fy2[ill]);

  *fwhp = xhpu - xhpl;
  return (0);
}
