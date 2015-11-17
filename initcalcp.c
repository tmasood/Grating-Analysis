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

// initcalcp.c:  routines to initialize panel potentials

#include "gtoothpnl.h"

extern "C"
{
  void Jacobi(int, double, double, double*, double*);
}

void getlogweights(int, double*, double*);

double *xleg, *xlog, *wlog, *wleg; // Gauss abscissas/weights
int ordleg, ordlog;                // Gauss quadrature orders


void initcalcp(int qorder)
{
  ordleg = qorder;
  ordlog = qorder;

  // abscissas, weights of Gauss-Legendre, scaled for [0,1]
  xleg = new double[ordleg];
  wleg = new double[ordleg];
  Jacobi(ordleg, 0.0, 0.0, xleg, wleg);


  /* abscissas, weights of quadrule for [0,1] with log singul. at zero */
  xlog = new double[ordlog];
  wlog = new double[ordlog];
  getlogweights(ordlog, xlog, wlog);

}
