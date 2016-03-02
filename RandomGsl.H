//////////////////////////////////////////////////////////////////////
// Copyright 2014-2016 Jeffrey Comer
//
// This file is part of DiffusionFusion.
//
// DiffusionFusion is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
//
// DiffusionFusion is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with DiffusionFusion. If not, see http://www.gnu.org/licenses/.
/////////////////////////////////////////////////////////////////////////
// Author: Jeff Comer <jeffcomer at gmail> 
// A wrapper for Gnu Scientific Library random number generators
#ifndef RANDOMGSL_H
#define RANDOMGSL_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

class Random {

private:
  gsl_rng* gslRando;

public:

  // default constructor
  Random() {
    gslRando = gsl_rng_alloc(gsl_rng_mt19937);
    init(0);
  }

  // constructor with seed
  Random(unsigned long seed) {
    gslRando = gsl_rng_alloc(gsl_rng_mt19937);
    init(seed);
  }

  ~Random() {
    gsl_rng_free(gslRando);
  }

  // reinitialize with seed
  void init(unsigned long seed) {
    gsl_rng_set(gslRando, seed);
  }

  // return a number uniformly distributed between 0 and 1
  double uniform() {
    return gsl_rng_uniform(gslRando);
  }

  long poisson(double lambda) {
     return gsl_ran_poisson(gslRando, lambda);
  }

  // return a number from a standard gaussian distribution
  double gaussian() {
    return gsl_ran_ugaussian(gslRando);
  }

  // return a random long
  long integer() {
    return gsl_rng_get(gslRando);
  }

  // Student's t-distribution.
  double student(double nu) {
    return gsl_ran_tdist(gslRando, nu);
  }

};

#endif
