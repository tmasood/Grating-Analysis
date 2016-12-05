//
//  Copyright 2015, 2016 Quantum Designs LLC, Taha Masood,
//  Johannes Tausch and Jerome Butler
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
#include <fstream>
#include <sstream>
#include <cstdlib>

#include "gtoothpnl.h"
#include "grating.h"

using namespace std;

// Assign indices to panels. Corresponding panels on Gamma_0 and Gamma_2
// get the same index. Panels on Gamma_3 and Gamma_4 have no defined idxU.
// For Gamma_0 panels the index of the matching Gamma_2 panel is written
// in the array match[].

const double IPEPS = 1e-12;

int indxpnls(grating *gratptr, int *match)
{
  gtpanel *pnlptr;
  gtpanel *pnl, *pnl2;

  int i;
  int k;
  int index, type, type2;

  double d1, d2, dist;
  double *p1x0, *p2x0, *p1x1, *p2x1;

  pnlptr = gratptr->gtpnlptr;
  index = 0;

  for (i=0, pnl=pnlptr; pnl!=NULL; pnl=pnl->nextptr, i++)
    {
      type = pnl->gettype();
      if (type == 0)
	{
	  for (k=0, pnl2=pnlptr; pnl2 != NULL; pnl2=pnl2->nextptr, k++)
	    {
	      type2 = pnl2->gettype();
	      if (type2 == 2)
		{
		  p1x0 = pnl->getx0();
		  p2x0 = pnl2->getx0();
		  d1  = (p1x0[0] - p2x0[0]);
		  dist = pow(d1,2.0);
		  p1x1 = pnl->getx1();
		  p2x1 = pnl2->getx1();
		  d2 = (p1x1[0] - p2x1[0]);
		  dist += pow(d2,2.0);
		  if ( dist < IPEPS )
		    {
		      // break from for (k=0 ...
		      break;
		    }
		}
	    }
	  if (pnl2 == NULL)
	    {
	      cout << "Panel x0 = " << p1x0[0] << " ";
	      cout << p1x0[1] << " x1 = " << p1x1[0] << " ";
	      cout << p1x1[1] << "has no corresponding panel "<< endl;
	      exit(1);
	    }
	  else
	    {
	      // assign same index to panel 0 and panel 2. Same periodic
	      // boundary condition
	      pnl->setidx(index);
	      pnl2->setidx(index);
	      index++;
	      match[i] = k;
	    }
	} // if type == 0

      type = pnl->gettype();
      if ((type == 1) || (type > 2))
	{
	  pnl->setidx(index);
	  index++;
	}
    } // for pnl - top level

  return 0; // end index panels
} 
