#include <complx.h>

complx cfm( complx r, double l )   /* retn the product of two dp complx nos. */
{
   complx dz;
   dz.re = l*r.re;
   dz.im = l*r.im;
   return( dz );
}
