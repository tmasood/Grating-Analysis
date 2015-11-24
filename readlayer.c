//
//  Copyright 2015  Taha Masood, Johannes Tausch and Jerome Butler
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
#include <fstream>
#include <sstream>

#include "structure.h"
#include "layer.h"

using namespace std;

int readlayer(string layerbuf, structure *epiptr)
{
  layer *layptr;
  layer *tmpptr;
  double indx; // layer refractive index
  double loss;  // layer loss (1/cm)
  double lt; // layer thickness (um)
  char layerinfo[MAXCHAR];

  istringstream ins(layerbuf);
  ins >> indx >> loss >> lt >> layerinfo;
  layptr = new layer(indx, loss, lt, layerinfo);

  tmpptr = epiptr->layerptr;
  if (tmpptr == NULL)
    {
      epiptr->layerptr = layptr;
    }
  else
    {
      while (tmpptr->nextptr != NULL)
	{
	  tmpptr = tmpptr->nextptr;
	}
      tmpptr->nextptr = layptr;
      layptr->prevptr = tmpptr;
    }

  return 0;
}
