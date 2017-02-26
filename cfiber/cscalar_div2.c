#include<stdio.h>
#include<math.h>
#include "f2c.h"

complex cscalar_div2(real m, complex a)
{
  complex temp;
  
  temp.r = m/a.r;
  if (a.i != 0.0)
    {
      temp.i = m/a.i;
    }
  else
    {
      temp.i = 0.0;
    }

  return (temp);

}
