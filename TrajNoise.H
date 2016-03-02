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
#ifndef TRAJNOISE_H
#define TRAJNOISE_H

#include <limits>
#include "useful.H"
#include "Field.H"
#include "PiecewiseLinear.H"
#include "PiecewiseCubic.H"
#include "Event.H"
#include "TrajCostComputer.H"

class TrajNoise : public TrajCostComputer {
private:
  const Piecewise1d* diffuse;
  const Piecewise1d* force;
  const Piecewise1d* sigma;
  int T, X, D, FB;

public:
  TrajNoise(const Field** ppf, IndexList fieldSel0, const Event* event0, int eventNum0, double kt, int leastLocal) :
    TrajCostComputer(ppf, fieldSel0, event0, eventNum0, kt, leastLocal, 4) {

    if (tcFieldNum != 3) {
      fprintf(stderr, "ERROR trajCost noise takes three fields: (0) diffusivity (1) force (2) sigma (3) \n");
      exit(-1);
    }

    diffuse = dynamic_cast<const Piecewise1d*>(fieldList[fieldSel.get(0)]);
    if (diffuse == NULL) {
       fprintf(stderr,"ERROR trajCost noise requires field type linear or cubic.\n");
       exit(-1);
    }
    force = dynamic_cast<const Piecewise1d*>(fieldList[fieldSel.get(1)]);
    if (force == NULL) {
       fprintf(stderr,"ERROR trajCost noise requires field type linear or cubic.\n");
       exit(-1);
    }
    sigma = dynamic_cast<const Piecewise1d*>(fieldList[fieldSel.get(2)]);
    if (force == NULL) {
       fprintf(stderr,"ERROR trajCost noise requires field type linear or cubic.\n");
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
      return String("sigma");
    default:
      return String("UNKNOWN");
    }
  }


  // The heart of the method: Eq 2 of Comer, Chipot, Gonzalez.
  // Plus the sigma of Türkcan et al. 2012
  inline double ccgCost(double dt, double dx, double frc, double dif, double gradDif, double sigma) {
    // Store the data that can be used to reconstruct gt.
    double gtNumer0 = dx - beta*dif*frc*dt - gradDif*dt;
    double gtVar0 = 2.0*(dif*dt + sigma*sigma);
    // Avoid race conditions.
    gtNumer = gtNumer0;
    gtVar = gtVar0;

    //printf("dt %g\n", dt);

    return 0.5*log(2.0*M_PI*gtVar0) + 0.5*(gtNumer0*gtNumer0)/gtVar0;
  }

  double eventCost(int e) {
    // Sum of forceBias and forceSys.
    double frc = event[e].var[FB] + force->computeVal(event[e].var[X]);
    double sig = sigma->computeVal(event[e].var[X]);
    // The diffusivity and its gradient.
    double dif, gradDif;
    diffuse->computeValGrad(event[e].var[X], dif, gradDif);
    // Give a huge result for negative diffusivities.
    if (dif <= 0.0) return std::numeric_limits<double>::max();

    // Add the cost of this event.
    double cost = ccgCost(event[e].del[T], event[e].del[D], frc, dif, gradDif, sig);

    // Check for NaN.
#ifdef FUSION_DEBUG
    //printf("e %d x %g dif %g %g %g\n", e, event[e].var[X], diffuse->computeVal(1.0),  diffuse->computeZeroOrder(1.0), diffuse->get(12));
    nanCheck(cost,'e',e,event[e].del[D],frc,gradDif,dif);
#endif
    return cost;
  }
};

#endif
