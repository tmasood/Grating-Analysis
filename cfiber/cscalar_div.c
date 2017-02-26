#include<stdio.h>
#include<math.h>
#include "f2c.h"

complex cscalar_div(complex a, real m)
{
  complex temp;
  
  temp.r = a.r/m;
  temp.i = a.i/m;

  return (temp);

}
