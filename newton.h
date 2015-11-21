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
#include <iostream>

#ifndef NEWTON_INCLUDE
#define NEWTON_INCLUDE

/* Structure of Newton search */
typedef struct{
  double rcond; /* Reciprocal Matrix Condition */
  complex<double> z;
} matnew;

#endif
