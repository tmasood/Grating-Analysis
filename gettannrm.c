//
//  Copyright 2015 Quantum Designs LLC, Taha Masood, Johannes Tausch
//  and Jerome Butler
//
//  Permission to use, copy, and distribute this software and its
//  documentation for any purpose with or without fee is hereby granted,
//  provided that the above copyright notice appear in all copies and
//  that both the copyright notice and this permission notice appear
//  in supporting documentation.
//
//  This software is provided "as is" without express or implied warranty
//  to the extent permitted by applicable law.
//

#include <iostream>
#include <cstdlib>
#include "gtoothdmn.h"
#include "grating.h"

// get the tangent and the normal vectors
// the tangent is from x0 to x1
// on \Gamma_0 and \Gamma_2 the normal points outward of the strip
// on \Gamma_1 the normal points into the domain on the left
// on \Gamma_3 the normal points to the left
// on \Gamma_4 the normal points to the right

void gettannrm(grating *gratptr)
{
  int i;
  int pnltype;

  gtpanel *pnlptr;
  
  double *x0, *x1;
  double xcoll0, xcoll1;
  double tan0, tan1;
  double nrm0, nrm1;
  double len;

  pnlptr = gratptr->gtrefpnlptr;

  for (i=0; pnlptr->nextptr != NULL; pnlptr = pnlptr->nextptr, i++)
    {
      x1 = pnlptr->getx1();
      x0 = pnlptr->getx0();
      xcoll0 = (x1[0] + x0[0])*0.5;
      xcoll1 = (x1[1] + x0[1])*0.5;
      pnlptr->setxcoll(xcoll0, xcoll1);
      len = pnlptr->getlength();
      tan0 = (x1[0] - x0[0])/len;
      tan1 = (x1[1] - x0[1])/len;

      pnlptr->settan(tan0, tan1);
      pnltype = pnlptr->gettype();

      switch (pnltype)
	{
	case 0:
	  nrm0 = 0.0;
	  nrm1 = -1.0;
	  break;

	case 2:
	  nrm0 = 0.0;
	  nrm1 = 1.0;
	  break;

	case 1:
	  nrm0 = -tan1;
	  nrm1 = tan0;
	  break;

	case 3:
	  nrm0 = -1.0;
	  nrm1 = 0.0;
	  break;
	
	case 4:
	  nrm0 = 1.0;
	  nrm1 = 0.0;
	  break;

	default:
	  cout << "Panel type not found -- gettannrm" << endl;
	  exit(1);
	}
      pnlptr->setnrm(nrm0, nrm1);
    }
}
