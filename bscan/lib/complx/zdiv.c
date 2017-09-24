#include <math.h>
#include <complx.h>

complx zdiv(complx a, complx b)
{
   double bm, cc, cd, zabs();
   complx zdiv_v;

   bm = 1.0e0/zabs( b );
   cc = b.re*bm;
   cd = b.im*bm;
   zdiv_v.re = (a.re*cc + a.im*cd)*bm;
   zdiv_v.im = (a.im*cc - a.re*cd)*bm;
   return( zdiv_v );
}
