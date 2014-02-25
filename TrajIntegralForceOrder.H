// Author: Jeff Comer <jeffcomer at gmail>
#ifndef TRAJINTEGRALFORCEORDER_H
#define TRAJINTEGRALFORCEORDER_H

#include "useful.H"
#include "Field.H"
#include "PiecewiseLinear.H"
#include "PiecewiseCubic.H"
#include "Event.H"
#include "TrajCostComputer.H"

class TrajIntegralForceOrder : public TrajCostComputer {
private:
  const Piecewise1d* diffuse;
  const Piecewise1d* force;
  int T, X, D, FB;
  // List of all leastLocalField nodes reached in each block and their neighbors.
  IndexList* blockNodeList;
  double* blockCurrCost;
  double* blockLastCost;
  int blockNum;

public:
  TrajIntegralForceOrder(const Field** ppf, IndexList fieldSel0, const Event* event0, int eventNum0, double kt, int leastLocal, int blockLen) :
    TrajCostComputer(ppf, fieldSel0, event0, eventNum0, kt, leastLocal, 4, blockLen), blockNodeList(NULL) {

    if (tcFieldNum != 2) {
      fprintf(stderr, "ERROR trajCost IntegralForceOrder takes two fields: (0) diffusivity (1) force\n");
      exit(-1);
    }

    diffuse = dynamic_cast<const Piecewise1d*>(fieldList[fieldSel.get(0)]);
    if (diffuse == NULL) {
      fprintf(stderr,"ERROR trajCost IntegralForceOrder requires field type linear or cubic.\n");
      exit(-1);
    }
    force = dynamic_cast<const Piecewise1d*>(fieldList[fieldSel.get(1)]);
    if (force == NULL) {
      fprintf(stderr,"ERROR trajCost IntegralForceOrder requires field type linear or cubic.\n");
      exit(-1);
    }

    T = 1;
    X = 0;
    D = 0;
    FB = 2;

    genBlockNodes();
    updateLocal();
    cloneLast();
  }

  ~TrajIntegralForceOrder() {
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
    default:
      return String("UNKNOWN");
    }
  }

double blockCost(int b) {
    const int e0 = b*blockLen;
    const int blockLen1 = blockLen-1;
    // Assume that the timestep is the same for the whole block.
    const double dt = event[e0].del[T];

    // The value and derivative of the integrand beta*D(x,t)*f(x,t) + grad D(x,t)
    double* integ = new double[blockLen];
    double dispSum = 0.0; 
    double varSum = 0.0;
    double meanSum = 0.0;
    for (int eb = 0; eb < blockLen; eb++) {
      int e = eb + e0;

      // Calculate the force and diffusivity at this time.
      double dif, gradDif;
      diffuse->computeValGrad(event[e].var[X], dif, gradDif);
      double frc = event[e].var[FB] + force->computeVal(event[e].var[X]);

      // Get each integrand value.
      integ[eb] = beta*dif*frc + gradDif;
      //printf("block %d eb %d detDx %g dx %g\n", b, eb, dt*integ[eb], event[e].del[D]);

      // Accumulate the displacement to avoid wrapping issues.
      // Note! We integrate only from 0 to blockLen-1.
      if (eb < blockLen1) {
	dispSum += event[e].del[D];
	varSum += dif;
	//meanSum += integ[eb];
      }
    }

    // Integrate by Simpson's rule.
    // Note! We integrate only from 0 to blockLen-1.
    meanSum = 0.5*(integ[0] + integ[blockLen1]); // First and last elements.
    // Trapezoid
    for (int eb = 1; eb < blockLen1; eb++) meanSum += integ[eb];
    meanSum *= dt;
    // Simpson's
    //for (int eb = 1; eb < blockLen1; eb+=2) meanSum += 4.0*integ[eb] + 2.0*integ[eb+1];
    //meanSum *= dt/3.0;
    
    // A crude zeroth-order approx.
    //meanSum *= dt;
    //meanSum = 0.5*dt*blockLen1*integ[0];

    delete[] integ;

    // These are the mean and variance sums.
    varSum *= 2.0*dt;
    double gtNumer0 = dispSum - meanSum;
    gtNumer = gtNumer0;
    gtVar = varSum;

    // Here is the integrated cost of the block
    double cost = 0.5*log(2.0*M_PI*varSum) + 0.5*(gtNumer0*gtNumer0)/varSum;

    // Check for NaN.
#ifdef FUSION_DEBUG
    nanCheck(cost,'b',b,dispSum,meanSum,meanSum,varSum);
#endif

    return cost;
  }

  void genBlockNodes() {
    // Sanity check.
    if (eventNum % blockLen != 0) {
      fprintf(stderr,"ERROR eventNum should be a multiple of blockLen!\n");
      exit(-1);
    }

    // Assign nodes to blocks.
    blockNum = eventNum/blockLen;
    blockNodeList = new IndexList[blockNum];
    blockCurrCost = new double[blockNum];
    blockLastCost = new double[blockNum];

#pragma omp parallel for
    for (int b = 0; b < blockNum; b++) {
      // Get the nodes for all events in this block.
      for (int eb = 0; eb < blockLen; eb++) {
	int e = eb + b*blockLen; // current event
	int near = leastLocalField->getNode(event[e].var);
	IndexList neigh = leastLocalField->neighbors(near);

	// We include both the nodes and their neighbors.
	for (int n = 0; n < neigh.length(); n++) {
	  int j = neigh.get(n);
	  int ind = blockNodeList[b].find(j);
	  // Add the node if it's not already in the list.
	  if (ind < 0) blockNodeList[b].add(j);
	} // Done with the neighbors of event eb of block b.
      } // Done with the events in block b.

      // Compact the memory used.
      blockNodeList[b].economize();
    } // End of block loop
  }

  // We return the number of blocks instead of the number of events.
  virtual int unitNum() const {
    return blockNum;
  }

  virtual double eventCost(int e) {
    // I used to always produce an error here, since we
    // only define costs for blocks of events.
    // However, it's simpler for some purposes to have
    // this give a block cost.
    // The function unitNum() ensures the caller iterates through
    // the blocks instead.

    if (e < 0 || e > blockNum) {
      fprintf(stderr,"ERROR eventCost called for TrajIntegralForceOrder with e greater than the number of blocks. TrajIntegralForceOrder gives the cost for a block of events, not for the events themselves.\n");
      exit(-1);
    }
    return blockCost(e);
  }


  ////////////////////////////////////////////////////////////////////
  // We overload calcCost, deltaCost, revert, cloneLast, and updateLocal.
  // The cost comes not from individual events, but from blocks.

  // Calculate the cost over all nodes.
  virtual double calcCost() {
    long double cost = 0.0;
#pragma omp parallel for reduction(+:cost)
    for (int b = 0; b < blockNum; b++) {
      cost += blockCost(b);
    }
    return cost;
  }

  virtual double deltaCost(const TrialMove& trialMove) {
    IndexList neigh = leastLocalField->neighbors(trialMove.node);

    // Find the blocks that are affected by a change in node.
    double dc = 0.0;
#pragma omp parallel for reduction(+:dc)
    for (int b = 0; b < blockNum; b++) {
      if ( blockNodeList[b].find(trialMove.node) >= 0) {
	// Block b is affected.
	// Save the last cost.
	blockLastCost[b] = blockCurrCost[b];
	// Calculate the new cost.
	double currCost = blockCost(b);
	// Set the current value.
	blockCurrCost[b] = currCost;
	// Add the difference.
	dc += currCost - blockLastCost[b]; 
      }
    }
    return dc;
  }

  virtual void revert(const TrialMove& trialMove) {
    // Find the blocks that are affected by a change in node.
#pragma omp parallel for
    for (int b = 0; b < blockNum; b++)
      if ( blockNodeList[b].find(trialMove.node) >= 0) {
	// Revert to the lastCost.
	blockCurrCost[b] = blockLastCost[b];
      }
  }

  virtual void cloneLast() {
    // Find the blocks that are affected by a change in node.
#pragma omp parallel for
    for (int b = 0; b < blockNum; b++)
	blockLastCost[b] = blockCurrCost[b];
  }

  // Instead of local costs, we have block costs.
  virtual double updateLocal() {
    long double cost = 0.0;
#pragma omp parallel for reduction(+:cost)
    for (int b = 0; b < blockNum; b++) {
      double bc = blockCost(b);
      blockCurrCost[b] = bc;
      cost += bc;
    }
    return cost;
  }
  
};

#endif
