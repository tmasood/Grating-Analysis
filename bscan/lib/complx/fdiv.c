#include <complx.h>

complx fdiv( a, b, c, d)
complx a, b, c, d;
{
   complx z;

   z = cdiv(tmul(a,b,c),d);
   return( z );
}   
