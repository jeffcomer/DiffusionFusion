/////////////////////////////////////////////////////////////////////////
// Author: Jeff Comer <jeffcomer at gmail> 
// Implementation of Equation 46 of
// Metzler and Klafter, Physics Reports 339 (2000) 1-77
// Using Gnu Scientific Library 

#ifndef FOXFUNC_H
#define FOXFUNC_H

#include <cmath>
#include <gsl/gsl_sf.h>

class FoxFunc {
private:
  int iterations;

public:
  // 40 iterations keeps errors below 5e-8 for alpha=0.1, larger alpha better
  FoxFunc() : iterations(40) {}
  FoxFunc(int iterations0) : iterations(iterations0) {}

  double stretchedGaussian(double x, double alpha, double sqrtKt, double fx) {
    double a1 = 1.0-alpha;
    double a2 = 2.0-alpha;
    double e = -0.5*a2*pow(0.5*alpha,alpha/a2)*pow(fx,1.0/(1.0-0.5*alpha));
    double res = sqrt(1.0/a2)*pow(2.0/alpha, a1/a2)*pow(fx,-a1/a2)*exp(e);
    return 0.5*res/(sqrtKt*sqrt(M_PI));
  }

  double w (double x, double alpha, double kAlpha, double t) {
    double sqrtKt = sqrt(kAlpha*pow(t,alpha));
    double fx = fabs(x)/sqrtKt;

    // Check if we can do the asymptotic instead.
    // The 20.0-12.66*alpha is based on empirical calculations
    // for getting a difference <1e-9 between the series and the asymptotic
    // difference of <1e-6, 12.108-7.55*alpha
    if (fx > 12.12 - 7.55*alpha) return stretchedGaussian(x, alpha, sqrtKt, fx);

    double sum = gsl_sf_gammainv(1.0-0.5*alpha);  
    int n = 1;
    double nFac = 1.0;
    double fxToN = 1.0;
    while (n < iterations) {
      nFac *= -n;
      fxToN *= fx;
      double term = gsl_sf_gammainv(1.0-0.5*alpha*(n+1))*fxToN/nFac;
      sum += term;
      n++;
    }

    return 0.5*sum/sqrtKt;
  }
};

#endif
