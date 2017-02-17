#include<stdio.h>
#include <stdlib.h>
#include<math.h>

int main(int argc, char *argv[])
{
  float a, b;
  float xn, x;
  float hv, f1, f2;
  float eo, edo, edd;
  float xo, xso, x1;
  float xlambda, x2, x3, x4;

  xlambda = atof(argv[1]);
  x = atof(argv[2]);

  a = (1 - x) * 5.14 + (x * 6.30);
  b = (1 - x) * 10.15 + (x * 9.40);
  hv = 1.24/xlambda;
  eo = 0.324 + (0.7 * x) + (0.4 * (x * x));
  edo = ((1.0 - x) * 0.38) + (x * 0.34);
  edd = edo + eo;
  xo = hv/eo;
  xso = hv/edd;
  f1 = (pow(xo,-2.0)) * (2 - (pow((1.0 + xo),0.5)) -
			 (pow(abs(1.0 - xo),0.5)));
  f2 = (pow(xso,-2.0)) * (2 - (pow((1 + xso),0.5)) - (pow(abs(1 - xso),0.5)));
  x1 = (pow((eo/edd),1.5))/2;
 printf (" xo = %f \n", xo);
  x2 = f1 + (x1 * f2);
  x3 = a * x2;
 printf ("x2 = %f, a = %f, x1 = %f \n",x2, a, x1);
  x4 = x3 + b;
 printf ("x4 = %f, x3 = %f, b = %f \n",x4, x3, b);
  xn = sqrt(x4);
  printf ("<th align=right> %f </th> \n", xn);
  return 0;
}
