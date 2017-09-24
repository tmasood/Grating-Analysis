#include <math.h>
#include <complx.h>

complx ccos( dz )   /* cosine of dp complx no. */
complx dz;
{
   complx dzc;

   dzc.re = cos(dz.re)*cosh(dz.im);
   dzc.im = -sin(dz.re)*sinh(dz.im);
   return( dzc );
}
