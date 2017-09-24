#include <complx.h>

Zmatrix Unit()
{
  Zmatrix v;
  (v.e11).re = 1.0;
  (v.e11).im = 0.0;
  (v.e21).re = 0.0;
  (v.e21).im = 0.0;
  v.e22 = v.e11;
  v.e12 = v.e21;
  return( v );
}  
