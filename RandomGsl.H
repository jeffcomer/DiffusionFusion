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
