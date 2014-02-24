///////////////////////////////////////////////////////////////////////
// Interface for a Smoluchowski equation solver
// Author: Jeff Comer <jeffcomer at gmail>
#ifndef SMOLUCHOWSKISOLVER_H
#define SMOLUCHOWSKISOLVER_H

#include "Piecewise1d.H"

class SmoluchowskiSolver {
public:
  virtual double solveProb(double x0, double x1, double tau, const Piecewise1d* diffuse, const Piecewise1d* forceSys, const Piecewise1d* forceBias, double& hopMax) const = 0;
  virtual void setParameters(double beta0, int steps0) = 0;
};

#endif
