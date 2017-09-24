#include <math.h>
#include <complx.h>

double zabs(complx z)
{
   double q, s, u, v, zabs_v;

   u = fabs( z.re );
   v = fabs( z.im );
   s = u + v;
   s *= 1.0e0;
   if( s == 0.0e0 ){
      zabs_v = 0.0e0;
      }
   else if( u > v ){
      q = v/u;
      zabs_v = u*sqrt( 1.e0 + q*q );
      }
   else{
      q = u/v;
      zabs_v = v*sqrt( 1.e0 + q*q );
      }
   return( zabs_v );
}
