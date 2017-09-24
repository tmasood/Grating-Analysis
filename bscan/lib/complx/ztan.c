#include <math.h>
#include <complx.h>

complx ztan(complx z)
{
   double t1, u, v, x, y;
   complx w;

   x =  z.re;
   y =  z.im;
   u = tan(x);
   v = tanh(y);
   t1 = 1.0 + u*u*v*v;
   w.re = u*(1.0 - v*v)/t1;
   w.im = v*(1.0 + u*u)/t1;
   return ( w );
} 
