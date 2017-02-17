#include <stdio.h>
#include "util.h"

extern double fcabs(double, double);

void clog(complex *r, complex *z)
{
  double zi, zr;
  r->i = atan2(zi = z->i, zr = z->r);
  r->r = log( fcabs(zr, zi) );
}
