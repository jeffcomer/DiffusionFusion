// Author: Jeff Comer <jeffcomer at gmail>
#ifndef TRAJREFLECT_H
#define TRAJREFLECT_H

#include "useful.H"
#include "Field.H"
#include "PiecewiseBicubic.H"
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
  TrajReflect(const Field** ppf, IndexList fieldSel0, const Event* event0, int eventNum0, double kt, int leastLocal) :
    TrajCostComputer(ppf, fieldSel0, event0, eventNum0, kt, leastLocal, 4) {

    if (tcFieldNum != 2) {
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
      return String("posX");
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
      return String("force");
    default:
      return String("UNKNOWN");
    }
  }

 // The heart of the method: Eq 2 of Comer , Chipot, Gonzalez.
  double reflectCost(double dt, double dx, double frc, double dif, double gradDif, double x0, double xa, double xb) {
    // Store the data that can be used to reconstruct gt.
    double mu = beta*dif*frc*dt + gradDif*dt;
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

  double eventCost(int e) {
    // Sum of forceBias and forceSys.
    double frc = event[e].var[FB] + force->computeVal(event[e].var[X]);
    // The diffusivity and its gradient.
    double dif = diffuse->computeVal(event[e].var[X]);
    double gradDif = diffuse->computeGrad(event[e].var[X]);

    // Add the cost of this event.
    return reflectCost(event[e].del[T], event[e].del[D], frc, dif, gradDif, event[e].var[D], bound0, bound1);
  }
};

#endif
