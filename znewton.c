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
#include <fstream>
#include <iostream>
#include <cstdlib>

#include "structure.h"
#include "newton.h"

using namespace std;

const int ITMAX = 50;
const double EPS = 1.0e-10;

matnew newtoncharmatrix(complex<double> *newzg, structure *epiptr);

complex<double> znewton(structure *epiptr)
{

  int i;
  int rootfound;

  double r1;
  double zguess_r;
  double zguess_i;

  matnew newtonmat;

  complex<double> zg;

  zguess_i = epiptr->getguessi();
  zguess_r = epiptr->getguessr();
  zg = complex<double>(zguess_r, zguess_i);

  newtonmat = newtoncharmatrix(&zg, epiptr);

  i = 1;

  while(i < ITMAX)
    {
      r1 = abs(newtonmat.z);
      if(r1 < EPS)
	{
	  rootfound = 1;
	  zg = zg + newtonmat.z;
	  return (zg);
	}
      else
	{
	  zg = zg + newtonmat.z;
	  newtonmat = newtoncharmatrix(&zg, epiptr);
	  i++;
	}
    }
  cout <<  "Failed to converge in " << ITMAX <<  " iterations." << endl;
  exit (1);
}
