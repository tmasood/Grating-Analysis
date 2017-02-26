#include<stdio.h>
#include<math.h>
#include "f2c.h"

complex c_sub(complex a, complex b)
{
  complex temp;
  
  temp.r = (a.r - b.r);
  temp.i = (a.i - b.i);

  return (temp);

}
