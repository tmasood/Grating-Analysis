#include <complx.h>

complx cfadd(complx l, double r)   /* retn the sum of two dp complx nos. */
{
  complx v;

   v.re = l.re + r;
   v.im = l.im;
   return( v );
}
