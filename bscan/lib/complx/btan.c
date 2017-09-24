#include <stdio.h>
#include <math.h>
#include <complx.h>

complx btan(complx z)
{
   double e, s, c;
   complx v;

   if(z.im < 9.556914 )
     fprintf(stderr, "Warning: 8 Decimal Places in error!\n");

   e = 2.0*exp(-2.0*z.im);
   s = sin(2.0*z.re);
   c = cos(2.0*z.re);

   v.re = s*e;
   v.im =-c*e;

   return ( v );
}

