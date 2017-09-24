#include <complx.h>

complx csub( complx l, complx r )   /* retn the difference of two dp complx nos. */
{
   l.re -= r.re;
   l.im -= r.im;
   return( l );
}
