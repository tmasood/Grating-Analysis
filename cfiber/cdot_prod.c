#include<stdio.h>
#include<math.h>
#include "f2c.h"

#define SIZE 4
complex cdot_prod(complex x, complex y)
{
  complex res;

  res.r = (x.r * y.r) + (x.i * y.i);
  res.i = 0.0;
  return (res);
}
