#include<stdio.h>
#include<math.h>
#include "f2c.h"
#include "cmplx.h"

int c_comparel(complex a, complex b)
{
  real x, y;
  
  x = c_abs(&a);
  y = c_abs(&b);
  if (x < y)
    {
      return (1);
    }
  else
    {
      return (0);
    }
}
