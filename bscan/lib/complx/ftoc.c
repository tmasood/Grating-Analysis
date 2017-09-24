#include <complx.h>

complx ftoc( r, i )   /* convert doubles to dp complx */
double r, i;
{
   complx dz;
   dz.re = r;   dz.im = i;
   return( dz );
}
