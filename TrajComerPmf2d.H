// Author: Jeff Comer <jeffcomer at gmail>
#ifndef TRAJCOMERPMF2D_H
#define TRAJCOMERPMF2D_H

#include "useful.H"
#include "Field.H"
#include "PiecewiseBicubic.H"
#include "Event.H"
#include "TrajCostComputer.H"

class TrajComerPmf2d : public TrajCostComputer {
private:
  const PiecewiseBicubic* diffuse;
  const PiecewiseBicubic* pmf;
  int T, X, Y, D, FB;

public:
  TrajComerPmf2d(const Field** ppf, IndexList fieldSel0, const Event* event0, int eventNum0, double kt, int leastLocal) :
    TrajCostComputer(ppf, fieldSel0, event0, eventNum0, kt, leastLocal, 5) {

    if (tcFieldNum != 2) {
      fprintf(stderr, "ERROR trajCost ccgPmf2d takes two fields: (0) diffusivity (1) pmf \n");
      exit(-1);
    }

    diffuse = dynamic_cast<const PiecewiseBicubic*>(fieldList[fieldSel.get(0)]);
    if (diffuse == NULL) {
       fprintf(stderr,"ERROR trajCost ccgPmf2d requires field type bicubic.\n");
       exit(-1);
    }
    pmf = dynamic_cast<const PiecewiseBicubic*>(fieldList[fieldSel.get(1)]);
    if (pmf == NULL) {
       fprintf(stderr,"ERROR trajCost ccgPmf2d requires field type bicubic.\n");
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
    for (int i = 0; i < 5; i++) printf("eventIndList[%d] = %d %s\n", i, eventIndList[i], eventVarName(i).cs());
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

    //printf("dt %g\n", dt);

    return 0.5*log(2.0*M_PI*gtVar0) + 0.5*(gtNumer0*gtNumer0)/gtVar0;
  }

  double eventCost(int e) {
    // Sum of forceBias and forceSys.
    double frc = event[e].var[FB] - pmf->computeGradX(event[e].var[X], event[e].var[Y]);
    // The diffusivity and its gradient.
    double dif = diffuse->computeVal(event[e].var[X], event[e].var[Y]);
    double gradDifX = diffuse->computeGradX(event[e].var[X], event[e].var[Y]);

    // Add the cost of this event.
    return ccgCost(event[e].del[T], event[e].del[D], frc, dif, gradDifX);
  }
};

#endif
