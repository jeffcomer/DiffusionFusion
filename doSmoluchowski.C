///////////////////////////////////////////////////////////////////////
// Author: Jeff Comer <jeffcomer at gmail>

#include <cstdio>
#include <cmath>

#include "useful.H"
#include "PiecewiseZero.H"
#include "PiecewiseCubic.H"

// 5-point numerical derivatives.
double firstDeriv5(double pm2, double pm1, double pc, double pp1, double pp2) {
  return (pm2 - 8.0*pm1 + 8.0*pp1 - pp2)/12.0;
}
double secondDeriv5(double pm2, double pm1, double pc, double pp1, double pp2) {
  return (-pm2 + 16.0*pm1 - 30.0*pc + 16.0*pp1 - pp2)/12.0;
}
// 3-point numerical derivatives.
double firstDeriv3(double pm1, double pc, double pp1) {
  return 0.5*(pp1-pm1);
}
double secondDeriv3(double pm1, double pc, double pp1) {
  return pm1 - 2.0*pc + pp1;
}

double slope(double dif, double gradDif, double frc, double gradFrc, double ddpi, double dpi, double pi) {
   return dif*ddpi + (gradDif - dif*frc)*dpi - (dif*gradFrc + gradDif*frc)*pi;
}


void finiteDif(const PiecewiseZero& p0, PiecewiseZero& p, double* dif, double* gradDif, double* frc, double* gradFrc, double dt, double dx) {
  double dx2 = dx*dx;

  // Get the 5-point stencil.
  int n = p0.length();

  // Absorbing boundary conditions at left.
  for (int i = 0; i < 2; i++) {
    double pm2 = 0.0;
    double pi = p0[i];
    double pp1 = p0[i+1];
    double pp2 = p0[i+2];
    double pm1;
    if (i == 0) pm1 = 0.0;
    else pm1 = p0[i-1];

    double dpi = firstDeriv5(pm2,pm1,pi,pp1,pp2)/dx;
    double ddpi = secondDeriv5(pm2,pm1,pi,pp1,pp2)/dx2;
    p.set(i, pi + slope(dif[i],gradDif[i],frc[i],gradFrc[i],ddpi,dpi,pi)*dt);
  }

  // Absorbing boundary conditions at right.
  for (int i = n-2; i < n; i++) {
    double pm2 = p0[i-2];
    double pm1 = p0[i-1];
    double pi = p0[i];
    double pp2 = 0.0;
    double pp1;
    if (i == n-1) pp1 = 0.0;
    else pp1 = p0[i+1];

    double dpi = firstDeriv5(pm2,pm1,pi,pp1,pp2)/dx;
    double ddpi = secondDeriv5(pm2,pm1,pi,pp1,pp2)/dx2;
    p.set(i, pi + slope(dif[i],gradDif[i],frc[i],gradFrc[i],ddpi,dpi,pi)*dt);
  }

  // The interior.
  for (int i = 2; i < n-2; i++) {
    double pm2 = p0[i-2];
    double pm1 = p0[i-1];
    double pi = p0[i];
    double pp1 = p0[i+1];
    double pp2 = p0[i+2];

    // 5-point numerical derivatives.
   double dpi = firstDeriv5(pm2,pm1,pi,pp1,pp2)/dx;
   double ddpi = secondDeriv5(pm2,pm1,pi,pp1,pp2)/dx2;
   p.set(i, pi + slope(dif[i],gradDif[i],frc[i],gradFrc[i],ddpi,dpi,pi)*dt);
  }
}


void finiteDif0(const PiecewiseZero& p0, PiecewiseZero& p, double* dif, double* gradDif, double* frc, double* gradFrc, double dt, double dx) {
  double dx2 = dx*dx;

  // Get the 5-point stencil.
  int n = p0.length();

  // Absorbing boundary conditions at left.
  {
    int i = 0;
    double pm1 = 0.0;
    double pi = p0[i];
    double pp1 = p0[i+1];

    double dpi = firstDeriv3(pm1,pi,pp1)/dx;
    double ddpi = secondDeriv3(pm1,pi,pp1)/dx2;
    p.set(i, pi + slope(dif[i],gradDif[i],frc[i],gradFrc[i],ddpi,dpi,pi)*dt);
  }

  // Absorbing boundary conditions at right.
  {
    int i = n-1;
    double pm1 = p0[i-1];
    double pi = p0[i];
    double pp1 = 0.0;

    double dpi = firstDeriv3(pm1,pi,pp1)/dx;
    double ddpi = secondDeriv3(pm1,pi,pp1)/dx2;
    p.set(i, pi + slope(dif[i],gradDif[i],frc[i],gradFrc[i],ddpi,dpi,pi)*dt);
  }

  // The interior.
  for (int i = 2; i < n-2; i++) {
    double pm1 = p0[i-1];
    double pi = p0[i];
    double pp1 = p0[i+1];

    double dpi = firstDeriv3(pm1,pi,pp1)/dx;
    double ddpi = secondDeriv3(pm1,pi,pp1)/dx2;
    p.set(i, pi + slope(dif[i],gradDif[i],frc[i],gradFrc[i],ddpi,dpi,pi)*dt);
  }
}

int main(int argc, char* argv[]) {
  if ( argc != 8 ) {
    printf("Usage: %s inFile diffuseFile forceFile dt time beta outFile\n", argv[0]);
    return 0;
  }
  const char* inFile = argv[1];
  const char* diffuseFile = argv[2];
  const char* forceFile = argv[3];
  double dt = strtod(argv[4], NULL);
  double tim = strtod(argv[5], NULL);
  double beta = strtod(argv[6], NULL);
  const char* outFile = argv[argc-1];

  // Load the diffusivity and force.
  PiecewiseCubic diffuse(diffuseFile, false);
  PiecewiseCubic force(forceFile, false);
  // Load the solution domain.
  PiecewiseZero dom0(inFile, false);
  PiecewiseZero dom(dom0);
  const int n = dom.length();
  double dx = dom.getDr();

  // Make pointers to the solution domain.
  PiecewiseZero* p0 = &dom0;
  PiecewiseZero* p = &dom;
  
  // Prepare to do the finite difference.
  int steps = int(ceil(tim/dt));

  // Make a list of the diffusivity and other values (maybe more efficient).
  double* arDif = new double[n];
  double* arGradDif = new double[n];
  double* arFrcBeta = new double[n];
  double* arGradFrcBeta = new double[n];
  for (int i = 0; i < n; i++) {
    double x = p->getNodePos(i);
    double dif, gradDif;
    diffuse.computeValGrad(x, dif, gradDif);
    double frc, gradFrc;
    force.computeValGrad(x, frc, gradFrc);

    arDif[i] = dif;
    arGradDif[i] = gradDif;
    arFrcBeta[i] = beta*frc;
    arGradFrcBeta[i] = beta*gradFrc;
  }

  for (int s = 0; s < steps; s++) {
    // Swap the pointers.
    PiecewiseZero* tmp = p0;
    p0 = p;
    p = tmp;

    // finite
    finiteDif0(*p0, *p, arDif, arGradDif, arFrcBeta, arGradFrcBeta, dt, dx);
  }

  p->write(outFile);

  delete[] arDif;
  delete[] arGradDif;
  delete[] arFrcBeta;
  delete[] arGradFrcBeta;
  return 0;
}

