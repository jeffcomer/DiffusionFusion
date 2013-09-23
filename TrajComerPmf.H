// Author: Jeff Comer <jeffcomer at gmail>
#ifndef TRAJCOMERPMF_H
#define TRAJCOMERPMF_H

#include "useful.H"
#include "Field.H"
#include "PiecewiseLinear.H"
#include "PiecewiseCubic.H"
#include "Event.H"
#include "TrajCostComputer.H"

class TrajComerPmf : public TrajCostComputer {
private:
  const Piecewise1d* diffuse;
  const Piecewise1d* pmf;
  int T, X, D, FB;

public:
  TrajComerPmf(const Field** ppf, IndexList fieldSel0, const Event* event0, int eventNum0, double kt, int leastLocal) :
    TrajCostComputer(ppf, fieldSel0, event0, eventNum0, kt, leastLocal, 4) {

    if (tcFieldNum != 2) {
      fprintf(stderr, "ERROR trajCost ccgPmf takes two fields: (0) diffusivity (1) pmf\n");
      exit(-1);
    }

    diffuse = dynamic_cast<const Piecewise1d*>(fieldList[fieldSel.get(0)]);
    if (diffuse == NULL) {
       fprintf(stderr,"ERROR trajCost ccgPmf requires field type linear or cubic.\n");
       exit(-1);
    }
    pmf = dynamic_cast<const Piecewise1d*>(fieldList[fieldSel.get(1)]);
    if (pmf == NULL) {
       fprintf(stderr,"ERROR trajCost ccgPmf requires field type linear or cubic.\n");
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
      return String("pmf");
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

    //printf("kT %g dt %g\n", 1.0/beta, dt);

    return 0.5*log(2.0*M_PI*gtVar0) + 0.5*(gtNumer0*gtNumer0)/gtVar0;
  }

  double eventCost(int e) {
    // Sum of forceBias and forceSys.
    double frc = event[e].var[FB] - pmf->computeGrad(event[e].var[X]);
    // The diffusivity and its gradient.
    double dif, gradDif;
    diffuse->computeValGrad(event[e].var[X], dif, gradDif);

    // Add the cost of this event.
    double cost = ccgCost(event[e].del[T], event[e].del[D], frc, dif, gradDif);
    if (cost != cost) {
      fprintf(stderr,"ERROR eventCost NaN\n");
      fprintf(stderr,"NAN event %d x %g dt %g dx %g frc %g dif %g\n",
	      e, event[e].var[X], event[e].del[T], event[e].del[D], frc, dif);
      exit(-1);
    }

    return cost;
  }
};

#endif
