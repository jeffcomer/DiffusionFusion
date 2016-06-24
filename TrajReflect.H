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
#ifndef TRAJREFLECT_H
#define TRAJREFLECT_H

#include "useful.H"
#include "Field.H"
#include "PiecewiseCubic.H"
#include "Event.H"
#include "TrajCostComputer.H"

//#include <gsl/gsl_sf_erf.h>

class TrajReflect : public TrajCostComputer {
private:
  const Piecewise1d* diffuse;
  const Piecewise1d* force;
  int T, X, D, FB;
  double bound0;
  double bound1;

public:
  TrajReflect(const TrajCostDesc& tcd) :
    TrajCostComputer(tcd, 4) {

    if (fieldSel.length() != 2) {
      fprintf(stderr, "ERROR trajCost reflect takes two fields: (0) diffusivity (1) force \n");
      exit(-1);
    }

    diffuse = dynamic_cast<const Piecewise1d*>(fieldList[fieldSel.get(0)]);
    if (diffuse == NULL) {
       fprintf(stderr,"ERROR trajCost reflect requires field type linear or cubic.\n");
       exit(-1);
    }
    force = dynamic_cast<const Piecewise1d*>(fieldList[fieldSel.get(1)]);
    if (force == NULL) {
       fprintf(stderr,"ERROR trajCost reflect requires field type linear or cubic.\n");
       exit(-1);
    }

    // Set the boundary based on the diffusivity field and the dimension.
    bound0 = diffuse->getR0();
    bound1 = diffuse->getR1();

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
    for (int i = 0; i < 4; i++) printf("eventIndList[%d] = %d %s\n", i, eventIndList[i], eventVarName(i).cs());
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
    default:
      return String("UNKNOWN");
    }
  }

  double reflectCost(double dt, double dx, double frc, double dif, double gradDif, double x0, double xa, double xb) {
    // Store the data that can be used to reconstruct gt.
    double mu = (beta*dif*frc + gradDif)*dt;
    double var = 2.0*dif*dt;
    double alpha = 0.5/var;
    
    double x1 = x0 + dx;
    double gtNumer0 = dx - mu;
    double plus =  x1 + x0 - mu;
    double numer1 = plus - 2.0*xa;
    double numer2 = plus + 2.0*xb;

    double sum = exp(-alpha*gtNumer0*gtNumer0);
    sum += exp(-alpha*numer1*numer1);
    sum += exp(-alpha*numer2*numer2);

    // Avoid race conditions.
    gtNumer = gtNumer0;
    gtVar = var;

    return 0.5*log(2.0*M_PI*var) - log(sum);
  }

  static void insertionSort(double* a, int len) {
    int i, j;
    double tmp;
    for (i = 1; i < len; i++) {
      j = i;
      while (j > 0 && a[j-1] > a[j]) {
	tmp = a[j];
	a[j] = a[j-1];
	a[j-1] = tmp;
	j--;
      }
    }
  }
   

  // Hopefully more numerically stable.
  double reflectCost1(double dt, double dx, double frc, double dif, double gradDif, double x0, double xa, double xb) {
    // Store the data that can be used to reconstruct gt.
    double mu = (beta*dif*frc + gradDif)*dt;
    double var = 2.0*dif*dt;
    double alpha = 0.5/var;
    
    double x1 = x0 + dx;
    double gtNumer0 = dx - mu;
    // The images move in the direction opposite the primary Gaussian.
    double plus =  x1 + x0 + mu;
    double numer1 = plus - 2.0*xa;
    double numer2 = plus - 2.0*xb;

    // Sort the exponents.
    double arg[3];
    arg[0] = alpha*gtNumer0*gtNumer0;
    arg[1] = alpha*numer1*numer1;
    arg[2] = alpha*numer2*numer2;
    insertionSort(arg, 3);
   
    // Now we should be able to handle this to avoid floating-point overflows or underflows.
    // arg0 is the smallest, so it dominates the probability.
    // The right term is log(1 + something_<=_one + something_<=_one).
    double cost = arg[0] - log(1.0 + exp(arg[0]-arg[1]) + exp(arg[0]-arg[2]));

    // Avoid race conditions.
    gtNumer = gtNumer0;
    gtVar = var;

    return 0.5*log(2.0*M_PI*var) + cost;
  }


  double eventCost(int e) {
    // Sum of forceBias and forceSys.
    double frc = event[e].var[FB] + force->computeVal(event[e].var[X]);
    // The diffusivity and its gradient.
    double dif = diffuse->computeVal(event[e].var[X]);
    double gradDif = diffuse->computeGrad(event[e].var[X]);

    // Add the cost of this event.
    return reflectCost1(event[e].del[T], event[e].del[D], frc, dif, gradDif, event[e].var[D], bound0, bound1);
  }
};

#endif
