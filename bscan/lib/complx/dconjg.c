#include <complx.h>

complx dconjg( complx z)
{
   complx w;
   w.re = z.re;
   w.im =-z.im;
   return( w );
}
