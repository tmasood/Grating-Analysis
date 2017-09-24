#include <math.h>
#include <complx.h>

complx zch(complx z)
{
   double ch, cn, sh, sn;
   complx zch_v;

   sh = sinh(z.re);
   ch = cosh(z.re);
   sn = sin(z.im);
   cn = cos(z.im);
   zch_v.re = ch*cn;
   zch_v.im = sh*sn;
   return( zch_v );
}
