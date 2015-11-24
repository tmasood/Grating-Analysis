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
#include "ssystem.h"

// default constructor
ssystem::ssystem()
{
  loopincr = 0.0;
}

// default destructor
ssystem::~ssystem()
{
}

void ssystem::setloopincr(double incr)
{
  loopincr = incr;
}

void ssystem::setnumsteps(int nsteps)
{
  numsteps = nsteps;
}

// counter propagating wave amplitude
void ssystem::setcpwamplitude(double amp)
{
  amplitude = amp;
}

// counter propagating wave phase
void ssystem::setcpwphase(double ph)
{
  phase = ph;
}
