#include<stdio.h>
#include<math.h>
#include "f2c.h"
#include "cmplx.h"

int c_comparel(complex a, complex b)
{
  if ((a.r < b.r) || ((a.r == b.r) && (a.i < b.i)))
    {
      return (1);
    }
  else
    {
      return (0);
    }  
}
