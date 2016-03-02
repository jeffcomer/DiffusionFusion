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
#ifndef TRAJLANGEVINEXP2D_H
#define TRAJLANGEVINEXP2D_H

#include "useful.H"
#include "Field.H"
#include "PiecewiseBicubic.H"
#include "Event.H"
#include "TrajCostComputer.H"

class TrajLangevinExp2d : public TrajCostComputer {
private:
  const PiecewiseBicubic* diffuse;
  const PiecewiseBicubic* pmf;
  const PiecewiseCubic* mass;
  int T, X, Y, D, FB;

public:
  TrajLangevinExp2d(const Field** ppf, IndexList fieldSel0, const Event* event0, int eventNum0, double kt, int leastLocal) :
    TrajCostComputer(ppf, fieldSel0, event0, eventNum0, kt, leastLocal, 5) {

    if (tcFieldNum != 3) {
      fprintf(stderr, "ERROR trajCost langevinExp2d takes three fields: (0) diffusivity (1) pmf (2) mass (a single entry)\n");
      exit(-1);
    }

    diffuse = dynamic_cast<const PiecewiseBicubic*>(fieldList[fieldSel.get(0)]);
    if (diffuse == NULL) {
       fprintf(stderr,"ERROR trajCost LangevinExp2d requires field type bicubic.\n");
       exit(-1);
    }
    pmf = dynamic_cast<const PiecewiseBicubic*>(fieldList[fieldSel.get(1)]);
    if (pmf == NULL) {
       fprintf(stderr,"ERROR trajCost langevinExp2d requires field type bicubic.\n");
       exit(-1);
    }
    mass = dynamic_cast<const PiecewiseCubic*>(fieldList[fieldSel.get(2)]);
    if (mass == NULL) {
       fprintf(stderr,"ERROR trajCost langevinExp2d requires field type cubic.\n");
       exit(-1);
    }

    T = 2;
    X = 0;
    Y = 1;
    D = 0;
    FB = 3;

    updateLocal();
    cloneLast();
  }

  // Event variables.
  void eventVarShortcuts() {
    T = eventIndList[0];
    X = eventIndList[1];
    Y = eventIndList[2];
    D = eventIndList[3];
    FB = eventIndList[4];
    for (int i = 0; i < 5; i++) printf("eventIndList[%d] = %d\n", i, eventIndList[i]);
  }
  // Event variables.
  String eventVarName(int ind) const {
    switch(ind) {
    case 0:
      return String("time");
    case 1:
      return String("posX");
    case 2:
      return String("posY");
    case 3:
      return String("displacement");
    case 4:
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
      return String("pmf");
    case 2:
      return String("mass");
    default:
      return String("UNKNOWN");
    }
  }


  // The heart of the method: Eq 2 of Comer, Chipot, Gonzalez.
  inline double lanExpCost(double dt, double dx, double dx0, double frc, double dif, double gradDif, double m) {
    double decay = exp(-dt/(m*beta*dif));
    double decay1 = 1.0 - decay;
    double gtNumer0 = dx - decay*dx0 - decay1*(beta*dif*frc + gradDif)*dt;
    double gtVar0 = (decay1*decay1)*(2.0*dif*dt);
    // Store the data that can be used to reconstruct gt.
    // Avoid race conditions.
    gtNumer = gtNumer0;
    gtVar = gtVar0;

    return 0.5*log(2.0*M_PI*gtVar0) + 0.5*(gtNumer0*gtNumer0)/gtVar0;
  }

  double eventCost(int e) {
    // Sum of forceBias and forceSys.
    double frc = event[e].var[FB] - pmf->computeGradX(event[e].var[X], event[e].var[Y]);
    // The diffusivity and its gradient.
    double dif = diffuse->computeVal(event[e].var[X], event[e].var[Y]);
    double gradDifX = diffuse->computeGradX(event[e].var[X], event[e].var[Y]);
    double m = mass->get(0);

    // Add the cost of this event.
    return lanExpCost(event[e].del[T], event[e].del[D], event[e].del0[D], frc, dif, gradDifX, m);
  }
};

#endif
