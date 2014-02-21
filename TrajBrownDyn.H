// Author: Jeff Comer <jeffcomer at gmail>
#ifndef TRAJBROWNDYN_H
#define TRAJBROWNDYN_H

#include "useful.H"
#include "Field.H"
#include "PiecewiseLinear.H"
#include "PiecewiseCubic.H"
#include "Event.H"
#include "TrajCostComputer.H"

struct BdParams {
  BdDescriptor(int stepsMin0, int stepsMax0, int hitsMin0) : stepsMin(stepsMin0), stepsMax(stepsMax0), hitsMin(hitsMin0) {
  }
public:
  int stepsMin;
  int stepsMax;
  int hitsMin;
};

class TrajBrownDyn : public TrajCostComputer {
private:
  const Piecewise1d* diffuse;
  const Piecewise1d* force;
  int T, X, D, FB;
  // Maximum node hop made by events in this node.
  int nodeHopMax*;
  Piecewise1d* bdBins;
  int bdStepsMin, bdStepsMax, bdHitsMin;

public:
  TrajBrownDyn(const Field** ppf, IndexList fieldSel0, const Event* event0, int eventNum0, double kt, int leastLocal, BdParams bdDesc0) :
    TrajCostComputer(ppf, fieldSel0, event0, eventNum0, kt, leastLocal, 4, blockLen), nodeHopMax(NULL), bdDesc(bdDesc0) {

    if (tcFieldNum != 3) {
      fprintf(stderr, "ERROR trajCost BrownDyn takes two fields: (0) diffusivity (1) force (2) brownianBins\n");
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
    Piecewise1d* bdBins0 = dynamic_cast<const Piecewise1d*>(fieldList[fieldSel.get(2)]);
    if (force == NULL) {
      fprintf(stderr,"ERROR trajCost BrownDyn requires field type linear or cubic.\n");
      exit(-1);
    }
    // The Brownian dynamics bins are a separate object with the same dimensions as the given field.
    bdBins = new Piecewise1d(*bdBins0);
    bdBins.zero();

    //std::numeric_limits<double>::infinity();

    T = 1;
    X = 0;
    D = 0;
    FB = 2;

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
      return String("force");
    case 2:
      return String("brownianBins");
    default:
      return String("UNKNOWN");
    }
  }


  void brownianDynamics() {
    
  }


  double eventCost(int e) {
    // We run Brownian dynamics on the current force and diffusivity fields
    // until our desired bin gets filled.
    for (int s = 0; s < bdStepsMax; s++) {
      
    }

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
