#include<stdio.h>
#include<math.h>
#include "f2c.h"
 
#define PI M_PI

#define SIZE 4

int assign_region(int nbound, real *xF, real *wrb, int *wregion)
{
  real bound;
  int i, kount;
  int region;
  
  kount = 200;
  region = 0;
  i = 0;
  bound = wrb[region + 1];
  while (i < kount)
    {
      if (xF[i] <= bound)
	{
	  wregion[i] = region;
	}
      else
	{
	  region++;
	  bound = wrb[region + 1];
	  wregion[i] = region;
	  if (region == nbound)
	    {
	      bound = 9999;
	    }
	}
      i++;
    }
  return (0);
}
