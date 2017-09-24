#include <complx.h>

complx cneg(complx dz)   /* dp complx unary minus operation */
{
   dz.re = -dz.re;   dz.im = -dz.im;
   return( dz );
}
