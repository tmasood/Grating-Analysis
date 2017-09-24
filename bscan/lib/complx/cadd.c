#include <complx.h>

complx cadd(complx l, complx r)   /* retn the sum of two dp complx nos. */
{
   l.re += r.re;
   l.im += r.im;
   return( l );
}
