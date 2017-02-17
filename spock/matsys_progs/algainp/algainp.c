#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[])
{
  float x, wvl, xn;
  float hv, eo, ed;

  int i;

  wvl = atof(argv[1]); 
  x = atof(argv[2]);

  hv = 1.24/wvl;
  eo = 3.39 + (0.62 * x);
  ed = 28.07 + (1.72 * x);
  xn = sqrt((eo * ed/((eo * eo) - (hv * hv))) + 1);
  printf ("<th align=right> %f </th> \n", xn);
  return (0);
}
