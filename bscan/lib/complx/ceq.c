#include <complx.h>

int ceq( l, r )   /* return TRUE if dp complx nos. equal */
complx l, r;
{
   if( l.re == r.re  &&  l.im == r.im )
      return( 1 );
   return( 0 );
}
