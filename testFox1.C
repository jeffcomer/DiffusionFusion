///////////////////////////////////////////////////////////////////////
// Author: Jeff Comer <jeffcomer at gmail>

#include <cstdio>
#include <cmath>
#include <gsl/gsl_sf.h>

double stretchedExp(double x, double alpha, double sqrtKt, double fx) {
  double a1 = 1.0-alpha;
  double a2 = 2.0-alpha;
  double e = -0.5*a2*pow(0.5*alpha,alpha/a2)*pow(fx,1.0/(1.0-0.5*alpha));
  double res = sqrt(1.0/a2)*pow(2.0/alpha, a1/a2)*pow(fx,-a1/a2)*exp(e);
  return 0.5*res/(sqrtKt*sqrt(M_PI));
}

double foxW(double x, double alpha, double kAlpha, double t, int iterations) {
  double sqrtKt = sqrt(kAlpha*pow(t,alpha));
  double fx = fabs(x)/sqrtKt;
  
  // Check if we can do the asymptotic instead.
  // The 20.0-12.66*alpha is based on empirical calculations
  // for getting a difference <1e-9 between the series and the asymptotic
  // difference of <1e-6, 12.108-7.55*alpha
  if (fx > 12.2 - 7.55*alpha) return stretchedExp(x, alpha, sqrtKt, fx);

  double sum = gsl_sf_gammainv(1.0-0.5*alpha);  
  int n = 1;
  double nFac = 1.0;
  double fxToN = 1.0;
  while (n < iterations) {
    double alt = (n%2==0)?1.0:-1.0;
    nFac *= n;
    fxToN *= fx;
    double term = alt*gsl_sf_gammainv(1.0-0.5*alpha*(n+1))*fxToN/nFac;
    sum += term;
    n++;
  }

  return 0.5*sum/sqrtKt;
}


int main(int argc, char* argv[]) {
  if ( argc != 7 ) {
    printf("Equation 45 & 46 of Metzler\n");
    printf("Usage: %s alpha t x0 x1 iterations outFile\n", argv[0]);
    return 0;
  }
  double alpha = strtod(argv[1], NULL);
  double t = strtod(argv[2], NULL);
  double x0 = strtod(argv[3], NULL);
  double x1 = strtod(argv[4], NULL);
  int iterations = atoi(argv[5]);
  const char* outFile = argv[argc-1];

  double dx = 0.01;
  double kAlpha = 1.0;

  FILE* out = fopen(outFile,"w");
  if (out == NULL) {
    fprintf(stderr,"ERROR Couldn't open file %s.\n",outFile);
    return -1;
  }

  for (double x = x0; x<= x1; x+=dx) {
    double y = foxW(x, alpha, kAlpha, t, iterations);
    fprintf(out, "%.15g %.15g\n", x, y);
  }
  fclose(out);

  return 0;
}
