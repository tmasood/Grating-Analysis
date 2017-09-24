#include <stdio.h>
#include <complx.h>
#include <math.h>

/*-------------------------------------------------------- 
 * This function evaluates the czeta function defined as
 *    \zeta(z) = tan(z) - j 
 *--------------------------------------------------------*/

complx czeta(complx z)
{
  double x, y, e, c, s;
  complx zeta;
  static complx J = {0.0,1.0};

  x = z.re;
  y = z.im;

  if( y < 9.55691 ) { /* 9.55691 from someplace */
    zeta = csub(ztan(z),J);
  }
  else{
    e = exp(-2.0*y);
    s = sin(x);
    c = cos(x);
    zeta.re = 4.0*c*s*e;
    zeta.im = -2.0*e;
  }
  return( zeta );
}
