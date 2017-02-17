#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[])
{
  float xx, qy, wvl;
  float xn;
  float a, b, c, x;
  float xlambda;

  wvl = atof(argv[1]); 
  xx = atof(argv[2]);
  qy = atof(argv[3]);

  xlambda = wvl * 1000.0;
  x = xx / 0.48;

  if((x < 0.3) && (wvl == 1.3))
    {
      xn = sqrt((1.0 - xx - qy)*14.6 + qy*13.2 + xx*10.06) - 0.17;
    }
  else if ((x  < 0.3) && (wvl  == 1.55))
    {
      xn = sqrt((1.0 - xx - qy)*14.6 + qy*13.2 + xx*10.06) - 0.31;
    }
  else if ((x < 0.3) && (wvl == 1.4))
    {
      xn = sqrt((1.0 - xx - qy)*14.6 + qy*13.2 + xx*10.06)- 0.226;
    }
  else if (x >= 0.3)
    {
      a = 9.689 - 1.012*x;
      b = 1.590 - 0.376*x;
      c = 1102.4 - 702.0*x + 330.4*(x*x);
      xn = sqrt(a + ((b*pow(xlambda,2))/(pow(xlambda,2) - pow(c,2))));
    }
  else
    {
      printf("please choose appropriate mole fractions and wavelength \n");
      exit(1);
    }
      
  printf ("<th align=right> %f </th> \n", xn);
  return (0);
}
 
