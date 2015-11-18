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
#include "layer.h"

// default constructor
layer::layer()
{
  index = 3.23;
  loss = 0.0;
  lt = 0.0;
  nlt = 0.0;
  kappa = complex<double>(0,0);
  hz = complex<double>(0,0);
  czh = complex<double>(0,0);
  szh = complex<double>(0,0);
  strcpy(comments,"Init");
  nextptr = NULL;
  prevptr = NULL;
}

layer::layer(double indx, double alpha, double th, char *layerinfo)
{
  index = indx;
  loss = alpha;
  lt = th;
  strcpy(comments,layerinfo);
  nextptr = NULL;
}

// default destructor
layer::~layer()
{
  delete [] nextptr;
}

void layer::setlayerkappa(complex<double> kap)
{
  kappa = kap;
}

void layer::setlayerkth(double nt)
{
   nlt = nt;
}

void layer::sethz(complex<double> ac)
{
   hz = ac;
}

void layer::setczh(complex<double> coszh)
{
   czh = coszh;
}

void layer::setszh(complex<double> sinzh)
{
   szh = sinzh;
}

double layer::getlayerindex()
{
  return index;
}

double layer::getlayerloss()
{
  return loss;
}

double layer::getlayerthickness()
{
  return lt;
}

double layer::getlayerkth()
{
  return nlt;
}

complex<double> layer::gethz()
{
   return hz;
}

complex<double> layer::getczh()
{
   return czh;
}

complex<double> layer::getszh()
{
   return szh;
}

complex<double> layer::getlayerkappa()
{
  return kappa;
}
