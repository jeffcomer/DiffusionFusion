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
  int T, X;
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
    
    leastLocalP1d = dynamic_cast<const Piecewise1d*>(leastLocalField);
    if (leastLocalP1d == NULL) {
      fprintf(stderr,"ERROR trajCost BdSmoluchowski requires leastLocalField type linear or cubic.\n");
      exit(-1);
    }

    // Set some physical parameters.
    solver->setParameters(beta,subintervalSteps);

    // For each node, we look for the neighboring nodes affected by events beginning in this node.
    deltaX = leastLocalP1d->getDr();
    nodeHopMax = new int[leastLocalNodes];
    for (int n = 0; n < leastLocalNodes; n++) nodeHopMax[n] = 3;
    updateHopMax();
    
    updateLocal();
    cloneLast();
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

  // Estimate the maximum hop distance for BD from each node.
  void updateHopMax() {
    /// In case we wanted a zero force.
    //PiecewiseCubic zero(*leastLocalP1d);
    //zero.zero();

    for (int n = 0; n < leastLocalNodes; n++) {
      int ne = local[n].events.length();

      // Find the biggest hop from each node.
      double hopMax = 0.0;
      for (int le = 0; le < ne; le++) {
	int e = local[n].events.get(le);
	double x0 = event[e].var[X];
	double x1 = x0 + event[e].del[X];
	const Piecewise1d* biasField = (event[e].bias < 0)?NULL:biasFieldList[event[e].bias];
	solver->solveProb(x0, x1, event[e].del[T], diffuse, force, biasField, hopMax);
      }

      // The +3 is because the change at i affects the
      // interpolation at i+2.
      int oldHopMax = nodeHopMax[n];
      nodeHopMax[n] = int(hopMax/deltaX)+3;
      //#ifdef FUSION_DEBUG
	if (abs(oldHopMax - nodeHopMax[n]) > 1) printf("nodeHopMax node %d changed from %d to %d\n", n, oldHopMax, nodeHopMax[n]);
      //#endif
    }
  }

  // Get the cost for an event by running Brownian Dynamics over the diffusivity and force fields.
  double eventCost(int e) {
    double hopMax = 0.0;
    double x0 = event[e].var[X];
    double x1 = x0 + event[e].del[X];
    const Piecewise1d* biasField = (event[e].bias < 0)?NULL:biasFieldList[event[e].bias];
    double prob = solver->solveProb(x0, x1, event[e].del[T], diffuse, force, biasField, hopMax);
    return -log(prob);
  }

  IndexList nodeRegion(int home) const {
    IndexList region;
    int n0 = home - nodeHopMax[home];
    int n1 = home + nodeHopMax[home];

    // BdSmoluchowskiSolver using diffuse, which could be
    // distinct from leastLocalField. TODO: YOU SHOULD CHECK FOR THIS!
    if (leastLocalField->getPeriodic()) {
      // Account for wrapping and make sure nodes are unique.
      for (int n = n0; n <= n1; n++) {
	int nw = leastLocalP1d->wrapNode(n);
	if (nw == n || nw < n0) region.add(n);
      }
    } else {
      if (n0 < 0) n0 = 0;
      if (n1 >= leastLocalNodes) n1 = leastLocalNodes-1;

      for (int n = n0; n <= n1; n++) region.add(n);
    }
    return region;
  }

  ////////////////////////////////////////////////////////////////////
  // We overload deltaCost, revert because we can no longer rely on just nearby neighbors.
  virtual double deltaCost(const TrialMove& trialMove) {
    IndexList region = nodeRegion(trialMove.node);
    return TrajCostComputer::deltaCost(trialMove, region);
  }

virtual void revert(const TrialMove& trialMove) {
    IndexList region = nodeRegion(trialMove.node);
    TrajCostComputer::revert(trialMove, region);
  }

  // We overload updateLocal() just so that we can update the maximum hops.
  virtual double updateLocal() {
    updateHopMax();
    return TrajCostComputer::updateLocal();
  }

};

#endif
