// Author: Jeff Comer <jeffcomer at gmail>
#ifndef TRAJBROWNDYN_H
#define TRAJBROWNDYN_H

#include "useful.H"
#include "Field.H"
#include "PiecewiseLinear.H"
#include "PiecewiseCubic.H"
#include "Event.H"
#include "TrajCostComputer.H"
#include "BrownDynDiffusionSolver.H"

class TrajBrownDyn : public TrajCostComputer {
private:
  const Piecewise1d* diffuse;
  const Piecewise1d* force;
  int T, X, D, FB;
  // Maximum node hop made by events in each node.
  int* nodeHopMax;
  DiffusionSolver* solver;

public:
  TrajBrownDyn(const Field** ppf, IndexList fieldSel0, const Event* event0, int eventNum0, double kt, int leastLocal, const Piecewise1d* biasFieldList, DiffusionSolver* solver0, int subintervalSteps) :
    TrajCostComputer(ppf, fieldSel0, event0, eventNum0, kt, leastLocal, 2), nodeHopMax(NULL), solver(solver0) {

    if (tcFieldNum != 2) {
      fprintf(stderr, "ERROR trajCost BrownDyn takes two fields: (0) diffusivity (1) force\n");
      exit(-1);
    }

    diffuse = dynamic_cast<const Piecewise1d*>(fieldList[fieldSel.get(0)]);
    if (diffuse == NULL) {
      fprintf(stderr,"ERROR trajCost BrownDyn requires field type linear or cubic.\n");
      exit(-1);
    }
    force = dynamic_cast<const Piecewise1d*>(fieldList[fieldSel.get(1)]);
    if (force == NULL) {
      fprintf(stderr,"ERROR trajCost BrownDyn requires field type linear or cubic.\n");
      exit(-1);
    }

    //std::numeric_limits<double>::infinity();

    T = 1;
    X = 0;
    
    nodeHopMax = new int[leastLocalField->length()];
    solver.setParameters(beta,subintervalSteps);

    updateLocal();
    cloneLast();
  }

  ~TrajBrownDyn() {
    if (blockNodeList != NULL) {
      delete[] blockNodeList;
      delete[] blockLastCost;
      delete[] blockCurrCost;
    }
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
    PiecewiseCubic zero(*leastLocalField);
    zero.zero();
    
    for (int e = 0; e < eventNum; e++) {
      hopMax = 0.0;
      solveProb(event[e].var[X], event[e].var[X]+event[e].del[X], event[e].del[T], diffuse, force, biasFieldList[event[e].bias], hopMax);
    }
  }
  
  double eventCost(int e) {
    
  }

  double eventCost(int e) {
    // Sum of forceBias and forceSys.
    double frc = event[e].var[FB] + force->computeVal(event[e].var[X]);
    // The diffusivity and its gradient.
    double dif, gradDif;
    diffuse->computeValGrad(event[e].var[X], dif, gradDif);

    // Add the cost of this event.
    double cost = ccgCost(event[e].del[T], event[e].del[D], frc, dif, gradDif);

    // Check for NaN.
#ifdef FUSION_DEBUG
    nanCheck(cost,'e',event[e].del[D],frc,gradDif,dif);
#endif
    return cost;
  }

  ////////////////////////////////////////////////////////////////////
  // We overload updateLocal, deltaCost, revert because we can no longer rely on just nearby neighbors.
  virtual double deltaCost(const TrialMove& trialMove) {
    
    
  }


  virtual double deltaCost(const TrialMove& trialMove) {
    IndexList neigh = leastLocalField->neighbors(trialMove.node);

    // Find the events that are affected by a change in node.
    double dc = 0.0;
#pragma omp parallel for reduction(+:dc)
    for (int b = 0; b < eventNum; b++) {
      if ( eventNodeList[b].find(trialMove.node) >= 0) {
	// Event b is affected.
	// Save the last cost.
	eventLastCost[b] = eventCurrCost[b];
	// Calculate the new cost.
	double currCost = eventCost(b);
	// Set the current value.
	eventCurrCost[b] = currCost;
	// Add the difference.
	dc += currCost - eventLastCost[b]; 
      }
    }
    return dc;
  }

  virtual void revert(const TrialMove& trialMove) {
    // Find the events that are affected by a change in node.
#pragma omp parallel for
    for (int b = 0; b < eventNum; b++)
      if ( eventNodeList[b].find(trialMove.node) >= 0) {
	// Revert to the lastCost.
	eventCurrCost[b] = eventLastCost[b];
      }
  }

  virtual void cloneLast() {
#pragma omp parallel for
    for (int b = 0; b < eventNum; b++)
	eventLastCost[b] = eventCurrCost[b];
  }

  // Instead of local costs, we have event costs.
  virtual double updateLocal() {
    long double cost = 0.0;
#pragma omp parallel for reduction(+:cost)
    for (int b = 0; b < eventNum; b++) {
      double bc = eventCost(b);
      eventCurrCost[b] = bc;
      cost += bc;
    }
    return cost;
  }
  
};

#endif
