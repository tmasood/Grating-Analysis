#include <stdio.h>
#include <math.h>
#include "util.h"
#include "layer.h"
#include "struct.h"

int ingaasp_adachi(struct LAYER *layerptr, struct STRUCT *structptr)
{
  float hv, eo, edo, xo, edd, xso;
  float a, b, fx, fxso, x1, x2, x3;

  hv = 1.24/structptr->wvl;
  a = 8.4 - (3.4 * layerptr->yperc);
  b = 6.6 + (3.4 * layerptr->yperc);
  eo = 1.35 - (0.72 * layerptr->yperc) +
    (0.12 * layerptr->yperc * layerptr->yperc);
  edd = 1.466 - (0.557 * layerptr->yperc) +
    (0.129 * layerptr->yperc * layerptr->yperc);
  xo = hv/eo;
  xso = hv/edd;                  
  fx = ((1/xo) * (1/xo)) * (2.0 - (sqrt(1.0 + xo))
			    - (sqrt(1 - xo)));
  fxso = ((1/xso) * (1/xso)) * (2.0 - (sqrt(1.0 + xso))
			      - (sqrt(1 - xso)));
  x1 = pow((eo/edd),1.5)/2;
  x2 = fx + (x1 * fxso);
  x3 = (a * x2) + b;
  layerptr->nreal = sqrt(x3);
  return (0);
}
 
