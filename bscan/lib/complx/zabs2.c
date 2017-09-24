#include <complx.h>
#include <math.h>

double zabs2(complx z)
{
  double u, v, zabs2_v;

  u = fabs( z.re );
  v = fabs( z.im );

  zabs2_v = v*v+u*u;

  return( zabs2_v );
}
