// Author: Jeff Comer <jeffcomer at gmail>
#ifndef BROWNIANDYNAMICS_H
#define BROWNIANDYNAMICS_H

#include "Piecewise1d.H"

class BrownianDynamics {
  static double step(double x0, const Piecewise1d* diffuse, const Piecewise1d* force, double , double beta) {
    double dif, gradDif;
    diffuse->computeValGrad(event[e].var[X], dif, gradDif);
    

  }
};

#endif
