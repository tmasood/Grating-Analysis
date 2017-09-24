#include <complx.h>

complx fcsub(double l, complx r) /* ret the sub of a complx from a double.*/
{
   r.re = l - r.re;
   return( r );
}
