//////////////////////////////////////////////////////////////////////
// Copyright 2014-2016 Jeffrey Comer
//
// This file is part of DiffusionFusion.
//
// DiffusionFusion is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
//
// DiffusionFusion is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with DiffusionFusion. If not, see http://www.gnu.org/licenses/.
// Author: Jeff Comer <jeffcomer at gmail>
#ifndef TRAJCOSTCOMPUTER_H
#define TRAJCOSTCOMPUTER_H

#include "TrajCostDesc.H"
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
  const IndexList fieldSel;
  const Event* event;
  const Field* leastLocalField;
  int leastLocalNodes;
  int trajVarMin;
  int* eventIndList;
  // Each node of each field has events associated with it.
  LocalCost* local;

  int eventStart;
  int eventEnd;
  int group;

public:
  // Store variables which can be used to reconstruct gt.
  // We do not calculate gt directly because we don't want to have to do an
  // extra sqrt() every step.
  double gtNumer; // displacement - mean
  double gtVar; // sigma^2

public:  
  TrajCostComputer(TrajCostDesc& tcd, int trajVarMin0)
    : beta(1.0/tcd.kbt), fieldList((const Field**)tcd.fieldList), fieldSel(tcd.fieldSel),
      event(tcd.event), leastLocalField(fieldList[tcd.leastLocal]), trajVarMin(trajVarMin0) {
    if (fieldSel.length() == 0) {
      fprintf(stderr, "ERROR TrajCostComputer must act on at least one field.\n");
      exit(-1);
    }
    leastLocalNodes = leastLocalField->length();

    group = tcd.group;
    if (group < 0) {
      // Just select all events.
      eventStart = 0;
      eventEnd = tcd.eventNum-1;
    } else {
      // Find the starting and ending events.
      eventStart = -1;
      eventEnd = -1;
      for (int e = 0; e < tcd.eventNum; e++) {
	if (eventStart < -1 && event[e].group == group) eventStart = e;
	if (eventStart != -1 && event[e].group != group) eventEnd = e-1;

	if (eventEnd != -1 && event[e].group == group) {
	  fprintf(stderr, "ERROR TrajCostComputer Events with the same group number must be contiguous.\n");
	  exit(-1);
	}
      }
      // This group reaches the end of the array.
      if (eventEnd == -1) eventEnd = tcd.eventNum-1;

      // Check that we found this group.
      if (eventStart == -1) {
	fprintf(stderr, "ERROR TrajCostComputer No events with -group %d\n", group);
	exit(-1);
      }
    }
        
    initLocal();
    eventIndList = new int[trajVarMin];
    for (int i = 0; i < trajVarMin; i++) eventIndList[i] = 0;
  }

  virtual ~TrajCostComputer() {
    // I'm not sure why we were deallocating event here, it's owned by DiffusionFusion.
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

  int getEventNum() const { return eventEnd - eventStart + 1; }
  int getTrajVarMin() const { return trajVarMin; }
  int eventVarIndex(int i) const {
    if (i >= 0 && i < trajVarMin) return eventIndList[i];
    return -1;
  }
  
  double getKt() const { return 1.0/beta; }
  int getNodes() const { return leastLocalNodes; }
  int getEventStart() const { return eventStart; }
  int getEventEnd() const { return eventEnd; }

  static bool nanCheck(double cost, int index, double disp, double mean, double mean1, double var) {
    if (cost != cost) {
      fprintf(stdout,"WARNING cost NaN\n");

      fprintf(stdout,"NAN event %d disp %g frc %g gradDif %g dif %g\n",index,disp,mean,mean1,var);
      return false;
    }
    return true;
  }

  // We'll allow the neighbors() to be overloaded in derived classes.
  // The default neighbors is derived from the number of nodes contributing
  // to the interpolant in the leastLocalField.
  virtual IndexList neighbors(int home) const {
    return leastLocalField->neighbors(home);
  }

  // Calculate the cost over all nodes.
  virtual double calcCost() {
    long double cost = 0.0;
#pragma omp parallel for reduction(+:cost)
    for (int e = eventStart; e <= eventEnd; e++) cost += eventCost(e);
    return cost;
  }

  // Calculate the change in cost using the LocalCost array.
  // Changes to a node in a field affect the events assigned
  // to that node, as well as those assigned to neighboring nodes.
  virtual double deltaCost(const TrialMove& trialMove) {
    IndexList neigh = neighbors(trialMove.node);
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
    IndexList neigh = neighbors(trialMove.node);
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
    for (int e = eventStart; e < eventEnd; e++) {
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
