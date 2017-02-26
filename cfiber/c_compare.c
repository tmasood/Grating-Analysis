#include<stdio.h>
#include<math.h>
#include "f2c.h"

int c_comparege(complex a, complex b)
{
  real x, y;
  int temp;
  
  x = c_abs(a);
  y = c_abs(b);
  if (x >= y)
    {
      return (1);
    }
  else
    {
      return (0);
    }
}
