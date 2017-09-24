#include <math.h>
#include <complx.h>

complx zexp(complx a)
{
   double zm;
   complx v;

   zm = exp( a.re );
   v.re = zm*cos(a.im);
   v.im = zm*sin(a.im);
   return( v );
}
