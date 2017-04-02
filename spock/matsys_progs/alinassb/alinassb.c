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
  a = ((1 - x) * y * 5.14) + ((1 - x) * (1-y) * 7.91) + 
    (x * y * 25.30) + (x * (1 - y) * 59.68);
  b = ((1 - x) * y * 10.15) + ((1 - x) * (1 - y) * 13.07) +
    (x * y * (-0.8)) + (x * (1 - y) * (-9.53));
  eo = 0.173 + (4.014 * x) - (1.418 * x * x);
  edo = ((1 - x) * y * 0.38) + ((1 - x) * (1 - y) * 0.81) + 
    (x * y * 0.28) + (x * (1 - y) * 0.65);
  edd = edo + eo;
  xo = hv/eo;
  xso = hv/edd;
  f1 = pow(xo,-2) * (2 - sqrt((1 + xo)) - sqrt((1 - xo)));
  f2 = pow(xso,-2) * (2 - sqrt((1 + xso)) - sqrt((1 - xso)));
  x1 = pow((eo/edd),1.5)/2;
  x2 = f1 + (x1 * f2);
  x3 = a * x2;
  x4 = x3 + b;
  xn = sqrt(x4);
  printf ("<th align=right> %f </th> \n", xn);
  return (0);
}

 
 