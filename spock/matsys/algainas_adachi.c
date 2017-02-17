#include <stdio.h>
#include <math.h>
#include "util.h"
#include "layer.h"
#include "struct.h"

int algainas_adachi(struct LAYER *layerptr, struct STRUCT *structptr)
{
  float hv, eo, edo, xo, edd, xso;
  float a, b, f1, f2, x1, x2, x3;
  float x4;

  hv = 1.24/structptr->wvl;
  eo = 0.75 + (1.548 * layerptr->xperc);
  edo = (layerptr->xperc * 0.28) + (layerptr->yperc * 0.34) +
    ((1 - layerptr->xperc - layerptr->xperc) * 0.38);
  xo = hv/eo;  
  edd = edo + eo;
  xso = hv/edd;                  
  a = (layerptr->xperc * 25.30) + (layerptr->xperc * 6.30) +
    ((1 - layerptr->xperc - layerptr->xperc)*5.14);
  b = (layerptr->xperc * (-0.80)) + (layerptr->xperc * 9.40) +
    ((1 - layerptr->xperc - layerptr->xperc) * 10.15);
  f1 = ((1/xo) * (1/xo)) * (2.0 - (sqrt(1.0 + xo))
			    - (sqrt(1 - xo)));
  f2 = ((1/xso) * (1/xso)) * (2.0 - (sqrt(1.0 + xso))
			      - (sqrt(1 - xso)));
  x1 = pow((eo/edd),1.5)/2;
  x2 = f1 + (x1 * f2);
  x3 = (a * x2) + b;
  layerptr->nreal = sqrt(x3);
  return (0);
}
 
