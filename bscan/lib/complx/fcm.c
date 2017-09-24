#include <complx.h>

complx fcm( double l, complx r )   /* retn the product of two dp complx nos. */
{
   complx dz;
   dz.re = l*r.re;
   dz.im = l*r.im;
   return( dz );
}
