// Author: Jeff Comer <jeffcomer at gmail>
#ifndef TRAJLANGEVINEXP_H
#define TRAJLANGEVINEXP_H

#include "useful.H"
#include "Field.H"
#include "PiecewiseCubic.H"
#include "Event.H"
#include "TrajCostComputer.H"

class TrajLangevinExp : public TrajCostComputer {
private:
  const Piecewise1d* diffuse;
  const Piecewise1d* pmf;
  const Piecewise1d* mass;
  int T, X, D, FB;

public:
  TrajLangevinExp(const Field** ppf, IndexList fieldSel0, const Event* event0, int eventNum0, double kt, int leastLocal) :
    TrajCostComputer(ppf, fieldSel0, event0, eventNum0, kt, leastLocal, 4) {

    if (tcFieldNum != 3) {
      fprintf(stderr, "ERROR trajCost langevinExp takes three fields: (0) diffusivity (1) pmf (2) mass (a single entry)\n");
      exit(-1);
    }

    diffuse = dynamic_cast<const Piecewise1d*>(fieldList[fieldSel.get(0)]);
    if (diffuse == NULL) {
       fprintf(stderr,"ERROR trajCost LangevinExp requires field type linear or cubic.\n");
       exit(-1);
    }
    pmf = dynamic_cast<const Piecewise1d*>(fieldList[fieldSel.get(1)]);
    if (pmf == NULL) {
       fprintf(stderr,"ERROR trajCost langevinExp requires field type linear or cubic.\n");
       exit(-1);
    }
    mass = dynamic_cast<const Piecewise1d*>(fieldList[fieldSel.get(2)]);
    if (mass == NULL) {
       fprintf(stderr,"ERROR trajCost langevinExp requires field type linear or cubic.\n");
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
    double frc = event[e].var[FB] - pmf->computeGrad(event[e].var[X]);
    // The diffusivity and its gradient.
    double dif, gradDif;
    diffuse->computeValGrad(event[e].var[X], dif, gradDif);
    double m = mass->get(0);

    // Add the cost of this event.
    return lanExpCost(event[e].del[T], event[e].del[D], event[e].del0[D], frc, dif, gradDif, m);
  }
};

#endif
