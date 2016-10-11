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

// calculate the RHS for 0-th order moments:
// the boundary condition on the left end first
// then the boundary condition on the right end

#include <iostream>
#include <complex>

#include "structure.h"
#include "layer.h"
#include "grating.h"

using namespace std;

// for a given point (x,z) in the interior domain return 
// in k the number of the domain.
// in alfa 1 or 1/2 if (x,z) is on the boundary or the interior
// The numbering starts at one, returns zero if (x,z) is not in 
// the interior.

extern double EPS;

// domainnr for sawtooth grating
int domainnrst(double x, double z, double *alf, structure *epiptr,
	     grating *gratptr)
{
  double thght;
  double period;
  double g1, g2;
  double slope;

  int gl;
  int layr;

  layer *layptr;

  *alf = 1.0;

  // get grating layer
  gl = gratptr->getgl();
  layptr = epiptr->layerptr;
  for (layr = 0; layr < gl; layr++)
    {
      layptr = layptr->nextptr;
    }

  thght = layptr->getlayerthickness();
  period = gratptr->getperiod();

  slope = 0.5*period/thght;

  if ( x<-EPS || x>thght+EPS || z<-EPS || z>period+EPS )
    {
      return 0;
    }

  g1 = z - slope*x;
  g2 = z - (1.0 - slope*x);

  if ( fabs(x) < EPS ) *alf = 0.5;
  if ( fabs(x - thght) < EPS ) *alf = 0.5;

  if ( fabs(z) < EPS ) *alf = 0.5;
  if ( fabs(g1) < EPS ) *alf = 0.5;
  if ( fabs(g2) < EPS ) *alf = 0.5;


  if ( g1 < 0.0) return 3;
  if (g2 > 0.0) return 1;

  return 2;
}
