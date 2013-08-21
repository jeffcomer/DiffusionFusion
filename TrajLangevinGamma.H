// Author: Jeff Comer <jeffcomer at gmail>
#ifndef TRAJLANGEVINGAMMA_H
#define TRAJLANGEVINGAMMA_H

#include "useful.H"
#include "Field.H"
#include "PiecewiseCubic.H"
#include "Event.H"
#include "TrajCostComputer.H"

class TrajLangevinGamma : public TrajCostComputer {
private:
  const PiecewiseCubic* friction;
  const PiecewiseCubic* pmf;
  const PiecewiseCubic* mass;
  int T, X, D, FB;

public:
  TrajLangevinGamma(const Field** ppf, IndexList fieldSel0, const Event* event0, int eventNum0, double kt, int leastLocal) :
    TrajCostComputer(ppf, fieldSel0, event0, eventNum0, kt, leastLocal, 4) {

    if (tcFieldNum != 3) {
      fprintf(stderr, "ERROR trajCost LangevinGamma takes three fields: (0) diffusivity (1) pmf (2) mass (a single entry)\n");
      exit(-1);
    }

    friction = dynamic_cast<const PiecewiseCubic*>(fieldList[fieldSel.get(0)]);
    if (friction == NULL) {
       fprintf(stderr,"ERROR trajCost LangevinGamma requires field type cubic.\n");
       exit(-1);
    }
    pmf = dynamic_cast<const PiecewiseCubic*>(fieldList[fieldSel.get(1)]);
    if (pmf == NULL) {
       fprintf(stderr,"ERROR trajCost LangevinGamma requires field type cubic.\n");
       exit(-1);
    }
    mass = dynamic_cast<const PiecewiseCubic*>(fieldList[fieldSel.get(2)]);
    if (mass == NULL) {
       fprintf(stderr,"ERROR trajCost LangevinGamma requires field type cubic.\n");
       exit(-1);
    }

    T = 0;
    X = 0;
    D = 0;
    FB = 0;

    updateLocal();
    cloneLast();
  }

  // Event variables.
  void eventVarShortcuts() {
    T = eventIndList[0];
    X = eventIndList[1];
    D = eventIndList[2];
    FB = eventIndList[3];
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
      return String("friction");
    case 1:
      return String("pmf");
    case 2:
      return String("mass");
    default:
      return String("UNKNOWN");
    }
  }


  // The heart of the method: Eq 2 of Comer, Chipot, Gonzalez.
  inline double bbkCost(double dt, double dx, double dx0, double frc, double gammaP, double m) {
    // Store the data that can be used to reconstruct gt.
    double dt2 = dt*dt;
    double denom = 1.0 + 0.5*gammaP*dt;
    double gtNumer0 = dx - ((1.0 - 0.5*gammaP*dt)*dx0 + dt2*frc/m)/denom;
    double gtVar0 = 2.0*gammaP*dt*dt2/(denom*denom*beta*m);
    // Avoid race conditions.
    gtNumer = gtNumer0;
    gtVar = gtVar0;

    return 0.5*log(2.0*M_PI*gtVar0) + 0.5*(gtNumer0*gtNumer0)/gtVar0;
  }

  double eventCost(int e) {
    // Sum of forceBias and forceSys.
    double frc = event[e].var[FB] - pmf->computeGrad(event[e].var[X]);
    double gammaP = friction->computeVal(event[e].var[X]);
    double m = mass->get(0);

    // Add the cost of this event.
    return bbkCost(event[e].del[T], event[e].del[D], event[e].del0[D], frc, gammaP, m);
  }
};

#endif
