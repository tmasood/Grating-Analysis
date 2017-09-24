#include <complx.h>

complx csq( complx l )   /* retn the square of a complx no. */
{
   complx dz;
   dz.re = l.re*l.re - l.im*l.im;
   dz.im = 2.0*l.re*l.im;
   return( dz );
}
