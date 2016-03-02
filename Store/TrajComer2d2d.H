// Author: Jeff Comer <jeffcomer at gmail>
#ifndef TRAJCOMER2D2D_H
#define TRAJCOMER2D2D_H

#include "useful.H"
#include "Field.H"
#include "PiecewiseBicubic.H"
#include "Event.H"
#include "TrajCostComputer.H"

class TrajComer2d2d : public TrajCostComputer {
private:
  const PiecewiseBicubic* diffuse;
  const PiecewiseBicubic* forceX;
  const PiecewiseBicubic* forceY;
  int T, X, Y, DX, DY, FBX, FBY;
  double invSqrt2;

public:
  TrajComer2d2d(const Field** ppf, IndexList fieldSel0, const Event* event0, int eventNum0, double kt, int leastLocal) :
    TrajCostComputer(ppf, fieldSel0, event0, eventNum0, kt, leastLocal, 7), invSqrt2(1.0/sqrt(2.0)) {

    if (tcFieldNum != 2) {
      fprintf(stderr, "ERROR trajCost ccg2d2d takes three fields: (0) diffusivity (1) forceX (2) forceY \n");
      exit(-1);
    }

    diffuse = dynamic_cast<const PiecewiseBicubic*>(fieldList[fieldSel.get(0)]);
    if (diffuse == NULL) {
       fprintf(stderr,"ERROR trajCost ccg2d2d requires field type bicubic.\n");
       exit(-1);
    }
    forceX = dynamic_cast<const PiecewiseBicubic*>(fieldList[fieldSel.get(1)]);
    if (forceX == NULL) {
       fprintf(stderr,"ERROR trajCost ccg2d2d requires field type bicubic.\n");
       exit(-1);
    }
    forceY = dynamic_cast<const PiecewiseBicubic*>(fieldList[fieldSel.get(1)]);
    if (forceY == NULL) {
       fprintf(stderr,"ERROR trajCost ccg2d2d requires field type bicubic.\n");
       exit(-1);
    }

    T = 2;
    X = 0;
    Y = 1;
    DX = 0;
    DY = 1;
    FBX = 3;
    FBY = 4;

    updateLocal();
    cloneLast();
  }

  // Event variables.
  void eventVarShortcuts() {
    T = eventIndList[0];
    X = eventIndList[1];
    Y = eventIndList[2];
    DX = eventIndList[3];
    DY = eventIndList[4];
    FBX = eventIndList[5];
    FBY = eventIndList[6];
    for (int i = 0; i < 7; i++) printf("eventIndList[%d] = %d\n", i, eventIndList[i]);
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
      return String("displacementX");
    case 4:
      return String("displacementY");
    case 5:
      return String("forceBiasX");
    case 6:
      return String("forceBiasY");
    default:
      return String("UNKNOWN");
    }
  }

  String fieldName(int ind) const {
    switch(ind) {
    case 0:
      return String("diffusivity");
    case 1:
      return String("forceX");
    case 2:
      return String("forceY");
    default:
      return String("UNKNOWN");
    }
  }

 // The heart of the method: Eq 2 of Comer, Chipot, Gonzalez.
  inline double ccgCost(double dt, double dx, double dy, double fx, double fy, double dif, double gradDifX, double gradDifY) {
    // Extra factor of two for cylindrical coordinates is probably approximate.
    double muX = dx - beta*dif*fx*dt - gradDifX*dt;
    double muY = dy - beta*dif*fy*dt - gradDifY*dt;
    // Store the data that can be used to reconstruct gt.
    gtNumer = (muX + muY)*invSqrt2;
    double gtVar0 = 2.0*dif*dt;
    // Avoid race conditions.
    gtVar = gtVar0;

    // Note that there's a factor of 2 on the first term.
    return log(2.0*M_PI*gtVar0) + 0.5*(muX*muX + muY*muY)/gtVar0;
  }

  double eventCost(int e) {
    // Sum of forceBias and forceSys.
    double fx = event[e].var[FBX] + forceX->computeVal(event[e].var[X], event[e].var[Y]);
    double fy = event[e].var[FBY] + forceY->computeVal(event[e].var[X], event[e].var[Y]);
    // The diffusivity and its gradient.
    double dif = diffuse->computeVal(event[e].var[X], event[e].var[Y]);
    double gradDifX = diffuse->computeGradX(event[e].var[X], event[e].var[Y]);
    double gradDifY = diffuse->computeGradY(event[e].var[X], event[e].var[Y]);

    // Add the cost of this event.
    return ccgCost(event[e].del[T], event[e].del[DX], event[e].del[DY], fx, fy, dif, gradDifX, gradDifY);
  }
};

#endif
