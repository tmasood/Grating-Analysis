#include<stdio.h>
#include<math.h>
#include "f2c.h"

complex c_prod(complex a, complex b)
{
  complex temp;
  
  temp.r = ((a.r*b.r) - (a.i*b.i));
  temp.i = ((a.r*b.i) + (a.i*b.r));

  return (temp);

}
