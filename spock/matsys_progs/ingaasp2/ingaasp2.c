#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[])
{
  float x, y, wvl;
  float xn;
  float hv, eo, edo, xo, edd, xso;
  float a, b, f1, f2, x1, x2, x3;
  float x4;

  wvl = atof(argv[1]); 
  x = atof(argv[2]);
  y = atof(argv[3]);

  hv = 1.24/wvl;
  a = 8.4 - (3.4 * y);
  b = 6.6 + (3.4 * y);
  eo = 1.35 - (0.72 * y) + (0.12 * y * y);
  edd = 1.466 - (0.557 * y) + (0.129 * y * y);
  xo = hv/eo;
  xso = hv/edd;                  
  if ((xo > 1.0) || (xso > 1.0))
   {
     xo = 1.0;
     xso = 1.0;
   } 
  f1 = ((1/xo) * (1/xo)) * (2.0 - (sqrt(1.0 + xo))
			    - (sqrt(1 - xo)));
  f2 = ((1/xso) * (1/xso)) * (2.0 - (sqrt(1.0 + xso))
			      - (sqrt(1 - xso)));
  x1 = pow((eo/edd),1.5)/2;
  x2 = f1 + (x1 * f2);
  x3 = a * x2;
  x4 = x3 + b;
  xn = sqrt(x4);
  printf ("<th align=right> %f </th> \n", xn);
  return (0);
}
 
