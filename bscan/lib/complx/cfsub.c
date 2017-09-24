#include <complx.h>

complx cfsub(complx l, double r) /* ret the sub of a double from complx. */
{
   l.re -= r;
   return( l );
}
