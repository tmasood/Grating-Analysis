#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[])
{
  float xx, qy, wvl;
  float xn;
  float hv, eo, edo, xo, edd, xso;
  float a, b, f1, f2, x1, x2, x3;
  float x4;

  wvl = atof(argv[1]); 
  xx = atof(argv[2]);
  qy = atof(argv[3]);

  hv = 1.24/wvl;
  eo = 0.75 + (1.548 * xx);
  edo = (xx * 0.28) + (qy * 0.34) + ((1 - xx - qy) * 0.38);
  xo = hv/eo;  
  edd = edo + eo;
  xso = hv/edd;                  
  a = (xx * 25.30) + (qy * 6.30) + ((1 - xx - qy)*5.14);
  b = (xx * (-0.80)) + (qy * 9.40) + ((1 - xx - qy) * 10.15);
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
 
