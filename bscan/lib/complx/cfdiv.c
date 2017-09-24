#include <complx.h>

complx cfdiv( complx r, double l )   /* retn the product of two dp complx nos. */
{
   complx dz;
   dz.re = r.re/l;
   dz.im = r.im/l;
   return( dz );
}
