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

#include "gtoothpnl.h"
#include "grating.h"
#include "structure.h"

using namespace std;

int defpnlst(structure *epiptr, grating *gratptr)
{
  layer *glayptr;
  gtpanel *pnlptr;

  int gl;
  int i;

  double thght;
  double period;

  glayptr = epiptr->layerptr;
  gl = gratptr->getgl();

  // point to the grating layer
  for (i = 0; i < gl; i++)
    {
       glayptr = glayptr->nextptr;
    }

  // get the grating thickness
  thght = glayptr->getlayerthickness();
  // get grating period
  period = gratptr->getperiod();
  // get the grating duty cycle

  // panel 0
  pnlptr = new gtpanel();
  pnlptr->setx0(0.0, 0.0);
  pnlptr->setx1(thght, 0.0);
  pnlptr->settype(0);
  gratptr->gtpnlptr = pnlptr;
  
  // panel 1
  pnlptr->nextptr = new gtpanel();
  pnlptr = pnlptr->nextptr;
  pnlptr->setx0(0.0, period);
  pnlptr->setx1(thght, period);
  pnlptr->settype(2);

  // panel 2
  pnlptr->nextptr = new gtpanel();
  pnlptr = pnlptr->nextptr;
  pnlptr->setx0(0.0, 0.0);
  pnlptr->setx1(thght, 0.5*period);
  pnlptr->settype(1);

  // panel 3
  pnlptr->nextptr = new gtpanel();
  pnlptr = pnlptr->nextptr;
  pnlptr->setx0(thght, 0.5*period);
  pnlptr->setx1(0.0, period);
  pnlptr->settype(1);

  // panel 4
  pnlptr->nextptr = new gtpanel();
  pnlptr = pnlptr->nextptr;
  pnlptr->setx0(0.0, 0.0);
  pnlptr->setx1(0.0, 0.5*period);
  pnlptr->settype(3);

  // panel 5
  pnlptr->nextptr = new gtpanel();
  pnlptr = pnlptr->nextptr;
  pnlptr->setx0(0.0, 0.5*period);
  pnlptr->setx1(0.0, period);
  pnlptr->settype(3);

  // panel 6
  pnlptr->nextptr = new gtpanel();
  pnlptr = pnlptr->nextptr;
  pnlptr->setx0(thght, 0.0);
  pnlptr->setx1(thght, 0.5*period);
  pnlptr->settype(4);

  // panel 7
  pnlptr->nextptr = new gtpanel();
  pnlptr = pnlptr->nextptr;
  pnlptr->setx0(thght, 0.5*period);
  pnlptr->setx1(thght, period);
  pnlptr->settype(4);

  // end of panel list
  pnlptr->nextptr = NULL;

  return 0;
}
