#include <complx.h>

complx fcadd(double l, complx r)   /* retn the sum of two dp complx nos. */
{
   r.re += l;
   return( r );
}
