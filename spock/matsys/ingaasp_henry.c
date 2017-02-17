#include<stdio.h>
#include<math.h>
#include "util.h"
#include "layer.h"
#include "struct.h"

int ingaasp_henry(struct LAYER *layerptr, struct STRUCT *structptr)
{
  float lambdapl;
  float h,a1,a2,ep,e1,e2,c;
  float epsilon, ee;
  float ec, term1, term2;

  lambdapl = layerptr->xperc;

  if (structptr->wvl <= lambdapl)
    {
      printf("<br><b>error: this material system not valid for \n");
      printf("lambda <= lambdapl\n</b><br>");
      exit(1);
    }

  if (structptr->wvl >= 1.8)
    {
      printf("<br><b>end wavelength should be less than \n");
      printf("1.8um \n</b><br>");
      exit(1);
    }

  /* Speed of light 'C', Plank's constant 'H'
     and electron charge 'EC' */
  c = 2.9979e8;
  h = 6.6261e-34;
  ec = 1.6022e-19;
  ep = (h*c)/(lambdapl*1e-6*ec);

  /* All energies are in Electron Volts */
  a1 = 13.3510 - (5.4554 * ep) + (1.2332 * pow(ep,2));
  a2 = 0.7140 - (0.3606 * ep);
          
  e1 = 2.5048;
  e2 = 0.1638;

  ee = (h * c)/(structptr->wvl * 1e-6 * ec);
  term1 = a1/(1 - (pow((ee/(ep + e1)),2)));
  term2 = a2/(1 - (pow((ee/(ep + e2)),2)));
  epsilon = 1 + term1 + term2;
  layerptr->nreal = sqrt(epsilon);

  return (0);
}
