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
#include "gtoothdmn.h"

// default constructor
gtdomain::gtdomain()
{
  npnls = 0;

  kap = complex<double>(0,0);
  indp = NULL;
  ornt = NULL;
  refindp = NULL;
  refornt = NULL;
  nextptr = NULL;
}

// default destructor
gtdomain::~gtdomain()
{
  delete [] indp;
  delete [] ornt;
  delete [] refindp;
  delete [] refornt;

  delete [] nextptr;
}

void gtdomain::setkap(complex<double> kappa)
{
  kap = kappa;
}

void gtdomain::setnpanels(int npanels)
{
  npnls = npanels;
}

void gtdomain::setnpts(int pts)
{
  npts = pts;
}

complex<double> gtdomain::getkap()
{
  return kap;
}

int gtdomain::getnpanels()
{
  return npnls;
}

int gtdomain::getnpts()
{
  return npts;
}
