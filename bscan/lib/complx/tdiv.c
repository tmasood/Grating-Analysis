#include <complx.h>

complx tdiv( l, c, r )   /* retn l*c/r complx nos. */
complx l, c, r;
{
   complx dz;

   dz = cmul(l,cdiv(c,r));
   return( dz );
}
