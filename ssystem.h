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
#include <iomanip>
#include "layer.h"

#ifndef SSYSTEM_INCLUDE
#define SSYSTEM_INCLUDE

class ssystem
{
  double loopincr;
  int numsteps;
  double amplitude;
  double phase;

 public:
  ssystem();
  ~ssystem();
  void setloopincr(double incr);
  void setnumsteps(int nsteps);
  void setcpwamplitude(double amp);
  void setcpwphase(double phase);
};

#endif
