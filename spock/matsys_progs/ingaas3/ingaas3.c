#include<stdio.h>
#include <stdlib.h>
#include<math.h>

/* Refractive index of Ga(1-x)In(x)As Prepared by Vapor-Phase Epitaxy */
/* by Toru Takagi */
/* Electrical communications lab., */
/* Nippon Telegraph and Telephone Public Corporation, */
/* Musahina-shi, Tokyo 180 */
/* October 1978 */

extern double pow(double, double);
extern double sqrt(double);

int main(int argc, char *argv[])
{
  double n, n2;
  double hv, x;
  double eo, ed, eg;
  double xlambda;

  xlambda = atof(argv[1]);
  x = atof(argv[2]);

  eg = 1.425 - (1.337 * x) + (0.27 * x * x); 
  if ((eg >= 0.5) && (eg <= 1.5))
    {
      hv = 1.24/xlambda;
      eo = 3.65 - (2.15 * x);
      ed = 36.1 - (19.9 * x);
      n2 = ((eo * ed) / (pow(eo,2.0) - pow(hv,2.0))) + 1.0;
      n = sqrt(n2);
      printf ("<th align=right> %f \t %f </th> \n", eg, n);
    }
  else
    {
      printf ("Model not valid for the specified wavelength and mole fraction \n");
    }
  return 0;
}
