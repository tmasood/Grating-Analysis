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

// default destructor
gtpanel::~gtpanel()
{
  delete [] nextptr;
}

void gtpanel::setx0(double x, double y)
{
  x0[0] = x;
  x0[1] = y;
}

void gtpanel::setx1(double x, double y)
{
  x1[0] = x;
  x1[1] = y;
}

void gtpanel::setxcoll(double x, double y)
{
  xcoll[0] = x;
  xcoll[1] = y;
}

void gtpanel::setlen(double length)
{
  len = length;
}

void gtpanel::settype(int typ)
{
  type = typ;
}

void gtpanel::setidx(int ind)
{
  idx = ind;
}

void gtpanel::setidxu(int indu)
{
  idxu = indu;
}

void gtpanel::setidxv(int indv)
{
  idxv = indv;
}

void gtpanel::settan(double t0, double t1)
{
  tan[0] = t0;
  tan[1] = t1;
}

void gtpanel::setnrm(double n0, double n1)
{
  nrm[0] = n0;
  nrm[1] = n1;
}

double *gtpanel::getx0()
{
  return x0;
}

double *gtpanel::getx1()
{
  return x1;
}

int gtpanel::gettype()
{
  return type;
}

int gtpanel::getidx()
{
  return idx;
}

int gtpanel::getidxu()
{
  return idxu;
}

int gtpanel::getidxv()
{
  return idxv;
}

double gtpanel::getlength()
{
  return len;
}

double *gtpanel::getxcoll()
{
  return xcoll;
}

double *gtpanel::gettan()
{
  return tan;
}

double *gtpanel::getnrm()
{
  return nrm;
}
