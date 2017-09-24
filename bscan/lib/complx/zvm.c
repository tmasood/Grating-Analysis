#include <complx.h>

zvector zvm(Zmatrix T, zvector v)
{
  zvector u;

  u.e1 = cadd(cmul(T.e11,v.e1),cmul(T.e12,v.e2));
  u.e2 = cadd(cmul(T.e21,v.e1),cmul(T.e22,v.e2));

  return( u );
}  
