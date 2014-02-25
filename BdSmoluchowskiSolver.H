///////////////////////////////////////////////////////////////////////
// Author: Jeff Comer <jeffcomer at gmail>
#ifndef BDSMOLUCHOWSKISOLVER_H
#define BDSMOLUCHOWSKISOLVER_H

#include "Piecewise1d.H"
#include "SmoluchowskiSolver.H"
#include "RandomGsl.H"

class BdSmoluchowskiSolver : public SmoluchowskiSolver {
public: 
  double binRad;
  int walksMin, walksMax, hitsMin;
  Random* rando;
  double posMin, posMax;

  double beta; // 1/(kT)
  int steps; // the number of subintervals of the propagation time
  
  BdSmoluchowskiSolver(double binSize0, int walksMin0, int walksMax0, int hitsMin0, Random* rando0) : binRad(binSize0*0.5), walksMin(walksMin0), walksMax(walksMax0), hitsMin(hitsMin0), rando(rando0), beta(-1.0) {}

  // This should be called by the TrajCostComputer, which has this information.
  virtual void setParameters(double beta0, int steps0) {
    beta = fabs(beta0);
    steps = abs(steps0);
  }

  // Get the probability of a hop from x0 to x1.
  // forceBias is allowed to be NULL, in which case the bias force is assumed to be zero.
  virtual double solveProb(double x0, double x1, double tau, const Piecewise1d* diffuse, const Piecewise1d* forceSys, const Piecewise1d* forceBias, double& hopMax) const {
    if (beta < 0.0) {
      fprintf(stderr,"ERROR BdSmoluchowskiSolver::solveProb You must call setParameters() first to set the physical parameters.\n");
      exit(-1);
    }

    // Run Brownian dynamics until hitsMin walkers land in the desired region or
    // we reach walksMax.
    double dt = tau/steps;
    int hits = 0;
    int w = 0;
    while( w <= walksMax && (w<walksMin || hits<hitsMin) ) {
      double xf = brownDyn(x0, dt, diffuse, forceSys, forceBias);
      double dx0 = fabs(xf-x0);
      double dx1 = fabs(xf-x1);
      // Periodicity is surmised from `diffuse'.
      if (diffuse->getPeriodic()) {
	dx0 = diffuse->wrapDiff(dx0);
	dx1 = diffuse->wrapDiff(dx1);
      }
      if (dx0 > hopMax) hopMax = dx0;

      // Did we land near the bin in question?
      if (dx1 < binRad) hits++;
      w++;
    }

    return hits/double(w);
  }


  // Run a set of Brownian dynamics steps.
  double brownDyn(double x0, double dt, const Piecewise1d* diffuse, const Piecewise1d* forceSys, const Piecewise1d* forceBias) const {
    double x = x0;

    for (int s = 0; s < steps; s++) {
      double frc = forceSys->computeVal(x);
      if (forceBias != NULL) frc += forceBias->computeVal(x);
      double dif, gradDif;
      diffuse->computeValGrad(x, dif, gradDif);
      x += (beta*frc*dif + gradDif)*dt + sqrt(2.0*dif*dt)*rando->gaussian();
    }
    
    // Periodicity is surmised from `diffuse'.
    if (diffuse->getPeriodic()) return diffuse->wrap(x);
    return x;
  }
    
};

#endif
