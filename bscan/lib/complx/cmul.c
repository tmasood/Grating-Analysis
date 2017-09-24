#include <complx.h>

complx cmul( complx l, complx r )   /* retn the product of two dp complx nos. */
{
   complx dz;
   dz.re = l.re*r.re - l.im*r.im;
   dz.im = l.re*r.im + l.im*r.re;
   return( dz );
}
