#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[])
{
  float x, wvl;
  float xn;
  float hv, eo, edo, xo, edd, xso;
  float a, b, f1, f2, x1, x2, x3;
  float x4;

  wvl = atof(argv[1]); 
  x = atof(argv[2]);

  hv = 1.24/wvl;
  a = ((1 - x) * 59.68) + (x * 25.30);
  b = ((1 - x) * (-9.53)) + (x * (-0.80));
  eo = 1.7 + (0.53 * x);
  edo = ((1 - x) * 0.65) + (x * 0.28);
  edd = (edo + eo);
  xo = hv/eo;
  xso = hv/edd;
  f1 = pow(xo,-2) * (2 - sqrt(1 + xo) - sqrt(1 - xo));
  f2 = pow(xso,-2) * (2 - sqrt(1 + xso) - sqrt(1 - xso));
  x1 = pow((eo/edd),1.5)/2;
  x2 = f1 + (x1 * f2);
  x3 = a * x2;
  x4 = x3 + b;
  xn = sqrt(x4);
  printf ("<th align=right> %f </th> \n", xn);
  return (0);
}
