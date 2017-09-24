#include <math.h>
#include <complx.h>

complx csin( dz )   /* sine of dp complx no. */
complx dz;
{
   complx dzs;

   dzs.re = sin( dz.re )*cosh( dz.im );
   dzs.im = cos( dz.re )*sinh( dz.im );
   return( dzs );
}
