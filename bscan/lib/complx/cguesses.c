#include <complx.h>

void cguesses(complx *z)
{
   *z = cadd(cmul(ftoc(0.99999,0.00001),*(z+2)),ftoc(0.00001,0.0));
   *(z+1) = csub(cmul(ftoc(1.00001,0.00001),*(z+2)),ftoc(0.00001,0.0));
}   
