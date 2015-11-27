/* Copyright 1985 Pierre Asselin.
/*
/* This package may be used and distributed freely, but may not be sold.
/* The present notice must be retained in all copies.
/* If incorporated in a commercial product,
/*  1) the package's origin and availability must be acknowledged
/*     prominently in the commercial product's documentation;
/*  2) the source code and documentation of the package must be
/*     made available to purchasers on request and at no extra cost.
/**/

/* Translation of module GaussKernel (implement) .
/**/

#include <stdio.h>		/* for error messages */
#include "GaussCore.h"

#define GaussEPS 1.0e-12
#define abs(x) ( ((x) < 0) ? -(x) : (x) )
#define sqr(x) ((x)*(x))

/*
/* Function to test a real value for exceeding -1,
/* and quit on an error if condition not met.
/**/
GaussCheck(value)
double value;
{
    if (value <= (-1.0)) {
	fprintf(stderr, "Gauss package: parameter out of range: %g\n", value);
	exit(1);
    }
}
