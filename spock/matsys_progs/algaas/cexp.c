#include <stdio.h>
#include "util.h"
 
void cexp(complex *r, complex *z)
{
  double expx, zi = z->i;
  
  expx = exp(z->r);
  r->r = expx * cos(zi);
  r->i = expx * sin(zi);
}
