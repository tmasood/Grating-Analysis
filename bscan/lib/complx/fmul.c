#include <complx.h>

complx fmul(a, b, c, d)
complx a, b, c, d;
{
   complx z;

   z = cmul(a,cmul(b,cmul(c,d)));
   return( z );
}
