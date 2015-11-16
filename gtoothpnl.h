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
// gtpanel - Grating tooth panel.
//
#include <iostream>

#define CSIZE 2

using namespace std;

#ifndef GTPANEL_INCLUDE
#define GTPANEL_INCLUDE

class gtpanel
{
  double x0[CSIZE];// start coordinate
  double x1[CSIZE]; // end coordinate
  double xcoll[CSIZE]; // x-collocation point
  double tan[CSIZE];
  double nrm[CSIZE];
  double len; // length of panel
  int type;   // panel type depending on which side it is facing
  int idx; // panel index
  int idxu; // index of u(xColl) in matrix (undefined for type 3&4)
  int idxv; // index of du/dn(xColl) in matrix


 public:
  gtpanel *nextptr;

  gtpanel();
  ~gtpanel();
  void setx0(double x, double y);
  void setx1(double x, double y);
  void setxcoll(double x, double y);
  void setlen(double length);
  void settype(int typ);
  void setidx(int ind);
  void setidxu(int indu);
  void setidxv(int indv);
  void settan(double t0, double t1);
  void setnrm(double n0, double n1);

  double *getx0();
  double *getx1();
  int gettype();
  int getidx();
  int getidxu();
  int getidxv();
  double getlength();
  double *getxcoll();
  double *gettan();
  double *getnrm();
};

#endif
