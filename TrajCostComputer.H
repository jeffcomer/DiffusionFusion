// Author: Jeff Comer <jeffcomer at gmail>
#ifndef TRAJCOSTCOMPUTER_H
#define TRAJCOSTCOMPUTER_H

#include "TrialMove.H"

struct LocalCost {
public:
  IndexList events;
  double lastCost;
  double currCost;
};

class TrajCostComputer {
protected:
  const double beta; // 1/kT
  // The fields that are involved in the trajCost.
  const Field** fieldList;
  const int tcFieldNum;
  const IndexList fieldSel;
  const int eventNum;
  const Event* event;
  const Field* leastLocalField;
  const int leastLocalNodes;
  int trajVarMin;
  int* eventIndList;
  int blockLen;
  // Each node of each field has events associated with it.
  LocalCost* local;

public:
  // Store variables which can be used to reconstruct gt.
  // We do not calculate gt directly because we don't want to have to do an
  // extra sqrt() every step.
  double gtNumer; // displacement - mean
  double gtVar; // sigma^2

public:  
  TrajCostComputer(const Field** ppf, IndexList fieldSel0, const Event* event0, int eventNum0, double kt, int leastLocalFieldInd, int trajVarMin0)
    : beta(1.0/kt), fieldList(ppf), tcFieldNum(fieldSel0.length()), fieldSel(fieldSel0), eventNum(eventNum0), event(event0), leastLocalField(ppf[leastLocalFieldInd]), leastLocalNodes(leastLocalField->length()), trajVarMin(trajVarMin0), blockLen(1) {
    if (tcFieldNum == 0) {
      fprintf(stderr, "ERROR TrajCostComputer must act on at least one field.\n");
      exit(-1);
    }
        
    initLocal();
    eventIndList = new int[trajVarMin];
    for (int i = 0; i < trajVarMin; i++) eventIndList[i] = 0;
  }

  TrajCostComputer(const Field** ppf, IndexList fieldSel0, const Event* event0, int eventNum0, double kt, int leastLocalFieldInd, int trajVarMin0, int blockLen0)
    : beta(1.0/kt), fieldList(ppf), tcFieldNum(fieldSel0.length()), fieldSel(fieldSel0), eventNum(eventNum0), event(event0), leastLocalField(ppf[leastLocalFieldInd]), leastLocalNodes(leastLocalField->length()), trajVarMin(trajVarMin0), blockLen(blockLen0) {
    if (tcFieldNum == 0) {
      fprintf(stderr, "ERROR TrajCostComputer must act on at least one field.\n");
      exit(-1);
    }

    if (blockLen < 1) {
      fprintf(stderr, "ERROR TrajCostComputer::blockLen must be at least 1.\n");
      exit(-1);
    }
    
    initLocal();
    eventIndList = new int[trajVarMin];
    for (int i = 0; i < trajVarMin; i++) eventIndList[i] = 0;
  }


  virtual ~TrajCostComputer() {
    delete[] event;
    delete[] local;
    delete[] eventIndList;
  }

  // Functions that must be implemented.
  // Note eventCost is not const because it sets gtNumer and gtVar.
  virtual double eventCost(int e) = 0;
  virtual String eventVarName(int ind) const = 0;
  virtual String fieldName(int ind) const = 0;
  virtual void eventVarShortcuts() = 0;

  void initCostVars(const String* costVarList, const IndexList& trajVarList) {
    for (int i = 0; i < trajVarList.length(); i++) {
      // Find the cost variable.
      int v;
      for (v = 0; v < trajVarMin; v++)
	if (eventVarName(v) == costVarList[i]) break;

      // Did we find a cost variable with the appropriate name?
      if (v == trajVarMin) {
	fprintf(stderr, "ERROR trajCost: Unrecognized trajCost variable name `%s'\n", costVarList[i].cs());
	fprintf(stderr, "Options are:");
	for (int j = 0; j < trajVarMin; j++) {
	  fprintf(stderr, " %s", eventVarName(j).cs());
	}
	fprintf(stderr, "\n");
	exit(-1);
      }
      
      eventIndList[v] = trajVarList.get(i);
      //printf("VARIABLE %s %d\n", costVarList[i].cs(), trajVarList.get(i));
    }

    eventVarShortcuts();
  }

  int getTrajVarMin() const { return trajVarMin; }
  int eventVarIndex(int i) const {
    if (i >= 0 && i < trajVarMin) return eventIndList[i];
    return -1;
  }
  
  double getKt() const { return 1.0/beta; }
  int getNodes() const { return leastLocalNodes; }
  // This is usually the number of events, but can
  // be overloaded to the number of blocks.
  virtual int unitNum() const { return eventNum; }

  static bool nanCheck(double cost, char unitType, int index, double disp, double mean, double mean1, double var) {
    if (cost != cost) {
      fprintf(stdout,"WARNING cost NaN\n");

      switch(unitType) {
      case 'e':
	fprintf(stdout,"NAN event %d disp %g frc %g gradDif %g dif %g\n",index,disp,mean,mean1,var);
	break;
      case 'b':
	fprintf(stdout,"NAN block %d disp %g mean %g var %g\n",index,disp,mean,var);
	break;

      }
      return false;
    }
    return true;
  }

  // Calculate the cost over all nodes.
  virtual double calcCost() {
    long double cost = 0.0;
#pragma omp parallel for reduction(+:cost)
    for (int e = 0; e < eventNum; e++) cost += eventCost(e);
    return cost;
  }

  // Calculate the change in cost using the LocalCost array.
  // Changes to a node in a field affect the events assigned
  // to that node, as well as those assigned to neighboring nodes.
  virtual double deltaCost(const TrialMove& trialMove) {
    IndexList neigh = leastLocalField->neighbors(trialMove.node);
    return deltaCost(trialMove,neigh);
  }
  virtual double deltaCost(const TrialMove& trialMove, const IndexList& neigh) {
      double dc = 0.0;
      // We got a >2 times speedup putting the parallel for here
      // instead of on the inner loop.
#pragma omp parallel for reduction(+:dc)
      for (int n = 0; n < neigh.length(); n++) {
	int j = neigh.get(n);
	// Save the current cost.
	local[j].lastCost = local[j].currCost;

	// Add the contributions of the events for the current cost.
	const int ne = local[j].events.length();
	double currCost = 0.0;
	for (int e = 0; e < ne; e++)
	  currCost += eventCost(local[j].events.get(e));

	// Set the current value.
	local[j].currCost = currCost;
	// Add the difference.
	dc += currCost - local[j].lastCost;
      }
    
      return dc;
  }

  inline int tcField(int id) const { return fieldSel.find(id); }

  // Revert a trial move.
  virtual void revert(const TrialMove& trialMove) {
    IndexList neigh = leastLocalField->neighbors(trialMove.node);
    revert(trialMove, neigh);
  }
  virtual void revert(const TrialMove& trialMove, const IndexList& neigh) {
#pragma omp parallel for
    for (int n = 0; n < neigh.length(); n++) {
      int j = neigh.get(n);
      // Revert to the lastCost.
      local[j].currCost = local[j].lastCost;
    }
  }

  // Set the local costs.
  virtual double updateLocal() {
    double cost = 0.0;
    
    for (int n = 0; n < leastLocalNodes; n++) {
      int ne = local[n].events.length();

      // Calculate the cost at each node.
      double currCost = 0.0;
#pragma omp parallel for reduction(+:currCost)
      for (int e = 0; e < ne; e++)
	currCost += eventCost(local[n].events.get(e));

      // Set the current value.
      local[n].currCost = currCost;
      // Accumulate the total cost.
      cost += currCost;
    }

    return cost;
  }

  // Make the last cost the current cost.
  virtual void cloneLast() {
     // Clone currCost to lastCost.
#pragma omp parallel for
    for (int n = 0; n < leastLocalNodes; n++)
      local[n].lastCost = local[n].currCost;
  }

  void printLocalCount() const {
      int count = 0;
      // Assign events to nodes.
      for (int n = 0; n < leastLocalNodes; n++)
	count += local[n].events.length();
      printf("count %d\n", count);
  }

private:
  // initLocal assumes that the ith coordinate trajectory variable
  // is aligned with the ith axis of the least local field.
  // Each event has a single home node.
  // Changes to a node in a field affect the events assigned
  // to that node, as well as those assigned to neighboring nodes.
  void initLocal() {
    // Set up the buffers for the cost of each node
    // of the least local field.
    local = new LocalCost[leastLocalNodes];

    // Assign events to nodes.
    for (int e = 0; e < eventNum; e++) {
      int near = leastLocalField->getNode(event[e].var);
      // Add this event to its near node.
      local[near].events.add(e);
    }

    // Economize memory after forming the event lists.
    for (int n = 0; n < leastLocalNodes; n++) local[n].events.economize();

    // You'll want to call updateLocal() after this.
    // However, eventCost() is implemented in the derived class,
    // so cannot be called (even indirectly) from the base class constructor.
  }

};
#endif
