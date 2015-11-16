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
#include "gtoothpnl.h"

// default constructor
gtpanel::gtpanel()
{
  x0[0] = 0.0;
  x0[1] = 0.0;
  x1[0] = 0.0;
  x1[1] = 0.0;
  xcoll[0] = 0.0;
  xcoll[1] = 0.0;
  tan[0] = 0.0;
  tan[1] = 0.0;
  nrm[0] = 0.0;
  nrm[1] = 0.0;
  len = 0.0;
  type = 0;
  idx = 0;
  idxu = 0;
  idxv = 0;
  
  nextptr = NULL;
}
