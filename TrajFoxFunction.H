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
// Author: Jeff Comer <jeffcomer at gmail>
#ifndef TRAJFRACTIONAL_H
#define TRAJFRACTIONAL_H

#include <limits>
#include "useful.H"
#include "Field.H"
#include "PiecewiseLinear.H"
#include "PiecewiseCubic.H"
#include "FoxFunc.H"
#include "Event.H"
#include "TrajCostComputer.H"

class TrajFoxFunction : public TrajCostComputer {
private:
  const Piecewise1d* diffuse;
  const Piecewise1d* force;
  const Piecewise1d* alpha; // Anomalous diffusion exponent
  int T, X, D, FB;
  FoxFunc fox;

public:
  TrajFoxFunction(const Field** ppf, IndexList fieldSel0, const Event* event0, int eventNum0, double kt, int leastLocal) :
    TrajCostComputer(ppf, fieldSel0, event0, eventNum0, kt, leastLocal, 4),
    fox(100) {

    if (tcFieldNum != 3) {
      fprintf(stderr, "ERROR trajCost fractional takes three fields: (0) diffusivity (1) force (2) anomalous_diffusion_exponent\n");
      exit(-1);
    }

    diffuse = dynamic_cast<const Piecewise1d*>(fieldList[fieldSel.get(0)]);
    if (diffuse == NULL) {
       fprintf(stderr,"ERROR trajCost fractional requires field type linear or cubic.\n");
       exit(-1);
    }
    force = dynamic_cast<const Piecewise1d*>(fieldList[fieldSel.get(1)]);
    if (force == NULL) {
       fprintf(stderr,"ERROR trajCost fractional requires field type linear or cubic.\n");
       exit(-1);
    }
    alpha = dynamic_cast<const Piecewise1d*>(fieldList[fieldSel.get(2)]);
    if (alpha == NULL) {
       fprintf(stderr,"ERROR trajCost fractional requires field type linear or cubic.\n");
       exit(-1);
    }

    T = 1;
    X = 0;
    D = 0;
    FB = 2;

    updateLocal();
    cloneLast();
  }

  // Event variables.
  void eventVarShortcuts() {
    T = eventIndList[0];
    X = eventIndList[1];
    D = eventIndList[2];
    FB = eventIndList[3];
    for (int i = 0; i < 4; i++) printf("eventIndList[%d] = %d\n", i, eventIndList[i]);
  }
  // Event variables.
  String eventVarName(int ind) const {
    switch(ind) {
    case 0:
      return String("time");
    case 1:
      return String("pos");
    case 2:
      return String("displacement");
    case 3:
      return String("forceBias");
    default:
      return String("UNKNOWN");
    }
  }

  String fieldName(int ind) const {
    switch(ind) {
    case 0:
      return String("diffusivity");
    case 1:
      return String("force");
    case 2:
      return String("alpha");
    default:
      return String("UNKNOWN");
    }
  }


  // The heart of the method: Eq 2 of Comer, Chipot, Gonzalez.
  inline double ccgCost(double dt, double dx, double frc, double dif, double gradDif) {
    // Store the data that can be used to reconstruct gt.
    double gtNumer0 = dx - beta*dif*frc*dt - gradDif*dt;
    double gtVar0 = 2.0*dif*dt;
    // Avoid race conditions.
    gtNumer = gtNumer0;
    gtVar = gtVar0;

    return 0.5*log(2.0*M_PI*gtVar0) + 0.5*(gtNumer0*gtNumer0)/gtVar0;
  }

  // The heart of the method: Equation 46 of Metzler and Klafter
  inline double fracCost(double dt, double dx, double frc, double dif, double gradDif, double alpha) {
    //printf("dt %g dx %g frc %g dif %g gradDif %g alpha %g\n", dt, dx, frc, dif, gradDif, alpha);
    double w = fox.w(dx, alpha, dif, dt);
    //if ( w != w || w<1e-300 ) printf("numerical problem with fox.w() %.15g\n", w);
    return -log(w);
  }

  double eventCost(int e) {
    // Sum of forceBias and forceSys.
    double frc = event[e].var[FB] + force->computeVal(event[e].var[X]);
    // The diffusivity and its gradient.
    double dif, gradDif;
    diffuse->computeValGrad(event[e].var[X], dif, gradDif);
    // Give a huge result for negative diffusivities.
    if (dif <= 0.0) return std::numeric_limits<double>::max();

    /// TESTING
    // Zero-out gradDif.
    gradDif = 0.0;

    // The anomalous diffusion exponent.
    double a = alpha->computeVal(event[e].var[X]);

    // Add the cost of this event.
    double cost = fracCost(event[e].del[T], event[e].del[D], frc, dif, gradDif, a);

    // Check for NaN.
#ifdef FUSION_DEBUG
    //printf("e %d x %g dif %g %g %g\n", e, event[e].var[X], diffuse->computeVal(1.0),  diffuse->computeZeroOrder(1.0), diffuse->get(12));
    nanCheck(cost,'e',e,event[e].del[D],frc,gradDif,dif);
#endif
    return cost;
  }
};

#endif
