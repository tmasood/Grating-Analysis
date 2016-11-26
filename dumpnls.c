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
#include "grating.h"

using namespace std;

// Refines the orignal panels.
// Returns an array of new panels

int dumpnls(grating *gratptr)
{
  gtpanel *pnlptr;
  int pnltype;
  int pnlidxu;
  int pnlidxv;
  double *pnlx0;
  double *pnlx1;
  double *pnlnrm;

  pnlptr = gratptr->gtrefpnlptr;
  cout << endl;
  cout << "Refine panels " << endl;
  while (pnlptr->nextptr != NULL)
    {
      cout << endl;
      pnltype = pnlptr->gettype();
      cout << "type = " << pnltype << endl;
      pnlx0 = pnlptr->getx0(); 
      cout << " x0 = " << pnlx0[0] << "    " << pnlx0[1];
      cout << endl;
      pnlx1 = pnlptr->getx1();
      cout << " x1 = " << pnlx1[0] << "   " << pnlx1[1];
      cout << endl;
      pnlnrm = pnlptr->getnrm();
      cout << " nrm = " << pnlnrm[0] << "   " << pnlnrm[1] << endl;
      pnlidxu = pnlptr->getidxu();
      pnlidxv = pnlptr->getidxv();
      cout << " idxu = " << pnlidxu << endl;
      cout << " idxv = " << pnlidxv << endl;
      cout << endl;

      pnlptr = pnlptr->nextptr;
    }

  return 0;
}
