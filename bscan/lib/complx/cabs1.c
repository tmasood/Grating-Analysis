#include <math.h>
#include <complx.h>

double cabs1(z)
complx z;
{
   double w;
   w = fabs(z.re) + fabs(z.im);
   return( w );
}   
