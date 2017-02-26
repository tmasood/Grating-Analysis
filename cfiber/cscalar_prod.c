#include<stdio.h>
#include<math.h>
#include "f2c.h"

complex cscalar_prod(complex a, real m)
{
  complex temp;
  
  temp.r = m*a.r;
  temp.i = m*a.i;

  return (temp);

}
