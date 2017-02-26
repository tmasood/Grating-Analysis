#include<stdio.h>
#include<math.h>
#include "f2c.h"
#include "cmplx.h"

int c_compareli(complex a, complex b)
{
  if ((a.i < b.i) || ((a.i == b.i) && (a.r < b.r)))
    {
      return (1);
    }
  else
    {
      return (0);
    }  
}
