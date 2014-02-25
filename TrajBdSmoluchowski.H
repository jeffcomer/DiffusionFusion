// Author: Jeff Comer <jeffcomer at gmail>
#ifndef TRAJBDSMOLUCHOWSKI_H
#define TRAJBDSMOLUCHOWSKI_H

#include "useful.H"
#include "Field.H"
#include "PiecewiseLinear.H"
#include "PiecewiseCubic.H"
#include "Event.H"
#include "TrajCostComputer.H"
#include "SmoluchowskiSolver.H"

class TrajBdSmoluchowski : public TrajCostComputer {
private:
  const Piecewise1d* diffuse;
  const Piecewise1d* force;
  const Piecewise1d** biasFieldList;
  int T, X;;
  // Maximum node hop made by events in each node.
  const Piecewise1d* leastLocalP1d;
  int* nodeHopMax;
  SmoluchowskiSolver* solver;
  double deltaX;

public:
  TrajBdSmoluchowski(const Field** ppf, IndexList fieldSel0, const Event* event0, int eventNum0, double kt, int leastLocal, const Piecewise1d** biasFieldList0, SmoluchowskiSolver* solver0, int subintervalSteps) :
    TrajCostComputer(ppf, fieldSel0, event0, eventNum0, kt, leastLocal, 2), biasFieldList(biasFieldList0), nodeHopMax(NULL), solver(solver0) {

    if (tcFieldNum != 2) {
      fprintf(stderr, "ERROR trajCost BdSmoluchowski takes two fields: (0) diffusivity (1) force\n");
      exit(-1);
    }

    diffuse = dynamic_cast<const Piecewise1d*>(fieldList[fieldSel.get(0)]);
    if (diffuse == NULL) {
      fprintf(stderr,"ERROR trajCost BdSmoluchowski requires field type linear or cubic.\n");
      exit(-1);
    }
    force = dynamic_cast<const Piecewise1d*>(fieldList[fieldSel.get(1)]);
    if (force == NULL) {
      fprintf(stderr,"ERROR trajCost BdSmoluchowski requires field type linear or cubic.\n");
      exit(-1);
    }

    //std::numeric_limits<double>::infinity();

    T = 1;
    X = 0;

    const Piecewise1d* leastLocalP1d = dynamic_cast<const Piecewise1d*>(leastLocalField);
    if (leastLocalP1d == NULL) {
      fprintf(stderr,"ERROR trajCost BdSmoluchowski requires leastLocalField type linear or cubic.\n");
      exit(-1);
    }
    deltaX = leastLocalP1d->getDr();

    nodeHopMax = new int[leastLocalNodes];
    solver->setParameters(beta,subintervalSteps);

    updateLocal();
    cloneLast();
    initMaxHop();
  }

  ~TrajBdSmoluchowski() {
    delete[] nodeHopMax;
  }

  // Event variables.
  void eventVarShortcuts() {
    T = eventIndList[0];
    X = eventIndList[1];
    for (int i = 0; i < 2; i++) printf("eventIndList[%d] = %d\n", i, eventIndList[i]);
  }
  // Event variables.
  String eventVarName(int ind) const {
    switch(ind) {
    case 0:
      return String("time");
    case 1:
      return String("pos");
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

  void initMaxHop() {
    PiecewiseCubic zero(*leastLocalP1d);
    zero.zero();

    for (int n = 0; n < leastLocalNodes; n++) {
      int ne = local[n].events.length();

      // Find the biggest hop from each node.
      double hopMax = 0.0;
      for (int le = 0; le < ne; le++) {
	int e = local[n].events.get(le);
	double x0 = event[e].var[X];
	double x1 = x0 + event[e].del[X];
	double tau = event[e].del[T];
	solver->solveProb(x0, x1, tau, diffuse, force, biasFieldList[event[e].bias], hopMax);
	solver->solveProb(x0, x1, tau, diffuse, &zero, biasFieldList[event[e].bias], hopMax);
      }

      nodeHopMax[n] = int(hopMax/deltaX)+1;
      printf("maxHop %d %g %d\n", n, hopMax, nodeHopMax[n]);
    }
  }

  void updateMaxHop() {
    for (int n = 0; n < leastLocalNodes; n++) {
      int ne = local[n].events.length();

      // Find the biggest hop from each node.
      double hopMax = 0.0;
      for (int le = 0; le < ne; le++) {
	int e = local[n].events.get(le);
	double x0 = event[e].var[X];
	double x1 = x0 + event[e].del[X];
	double tau = event[e].del[T];
	solver->solveProb(x0, x1, tau, diffuse, force, biasFieldList[event[e].bias], hopMax);
      }

      nodeHopMax[n] = int(hopMax/leastLocalP1d->getDr())+1;
    }
  }
  
  double eventCost(int e) {
    double hopMax = 0.0;
    return solver->solveProb(event[e].var[X], event[e].var[X]+event[e].del[X], event[e].del[T], diffuse, force, biasFieldList[event[e].bias], hopMax);
  }


  ////////////////////////////////////////////////////////////////////
  // We overload deltaCost, revert because we can no longer rely on just nearby neighbors.
  virtual double deltaCost(const TrialMove& trialMove) {
    IndexList region;
    int n0 = trialMove.node - nodeHopMax[trialMove.node];
    int n1 = trialMove.node + nodeHopMax[trialMove.node];
    for (int n = n0; n <= n1; n++) region.add(n);
    return TrajCostComputer::deltaCost(trialMove, region);
  }

virtual void revert(const TrialMove& trialMove) {
    IndexList region;
    int n0 = trialMove.node - nodeHopMax[trialMove.node];
    int n1 = trialMove.node + nodeHopMax[trialMove.node];
    for (int n = n0; n <= n1; n++) region.add(n);
    TrajCostComputer::revert(trialMove, region);
  }

  // We overload updateLocal() just so that we can update the maxHops.
  virtual double updateLocal() {
    updateMaxHop();
    return TrajCostComputer::updateLocal();
  }

};

#endif
