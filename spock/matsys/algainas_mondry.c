#include <stdio.h>
#include <math.h>
#include "util.h"
#include "layer.h"
#include "struct.h"

int algainas_mondry(struct LAYER *layerptr, struct STRUCT *structptr)
{
  float a, b, c, x;
  float xlambda;

  xlambda = structptr->wvl * 1000.0;
  x = layerptr->xperc / 0.48;

  if((x < 0.3) && (structptr->wvl == 1.3))
    {
      layerptr->nreal = 
	sqrt((1.0 - layerptr->xperc - layerptr->yperc)*14.6 +
	     layerptr->yperc*13.2 + layerptr->xperc*10.06) - 0.17;
    }
  else if ((x  < 0.3) && (structptr->wvl  == 1.55))
    {
      layerptr->nreal =
	sqrt((1.0 - layerptr->xperc - layerptr->yperc)*14.6 +
	     layerptr->yperc*13.2 + layerptr->xperc*10.06) - 0.31;
    }
  else if ((x < 0.3) && (structptr->wvl == 1.4))
    {
      layerptr->nreal =
	sqrt((1.0 - layerptr->xperc - layerptr->yperc)*14.6 +
	     layerptr->yperc*13.2 + layerptr->xperc*10.06) - 0.226;
    }
  else if (x >= 0.3)
    {
      a = 9.689 - 1.012*x;
      b = 1.590 - 0.376*x;
      c = 1102.4 - 702.0*x + 330.4*(x*x);
      layerptr->nreal =
	sqrt(a + ((b*pow(xlambda,2))/(pow(xlambda,2) - pow(c,2))));
    }
  else
    {
      printf("AlGaInAs-Mondry: Mole fraction not supported \n");
      exit(1);
    }
      
  return (0);
}
