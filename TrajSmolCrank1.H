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
#ifndef TRAJSMOLCRANK1_H
#define TRAJSMOLCRANK1_H

#include "useful.H"
#include "Field.H"
#include "PiecewiseLinear.H"
#include "PiecewiseCubic.H"
#include "PiecewiseZero.H"
#include "Event.H"
#include "TrajCostComputer.H"
#include "SolutionStruct.H"
#include "CrankNicholsonSolver1.H"

class TrajSmolCrank1 : public TrajCostComputer {
private:
  const Piecewise1d* diffuse;
  const Piecewise1d* force;
  const Piecewise1d* solnGrid;
  int T, X;
  double timestep, maxHop;
  
  SolutionStruct soln;
  PiecewiseCubic** solnProb;
  double delT;
  int steps;
  IndexList all;
  CrankNicholsonSolver1* solver;
  
public:
  TrajSmolCrank1(const Field** ppf, IndexList fieldSel0, const Event* event0, int eventNum0, double kt, int leastLocal, double timestep0, double maxHop0) :
    TrajCostComputer(ppf, fieldSel0, event0, eventNum0, kt, leastLocal, 2), timestep(timestep0), maxHop(maxHop0) {

    if (tcFieldNum != 3) {
      fprintf(stderr, "ERROR trajCost SmolCrank1 takes three fields: (0) diffusivity (1) force (2) solnGrid \n");
      exit(-1);
    }

    diffuse = dynamic_cast<const Piecewise1d*>(fieldList[fieldSel.get(0)]);
    if (diffuse == NULL) {
      fprintf(stderr,"ERROR trajCost SmolCrank1 requires field type linear or cubic.\n");
      exit(-1);
    }
    force = dynamic_cast<const Piecewise1d*>(fieldList[fieldSel.get(1)]);
    if (force == NULL) {
      fprintf(stderr,"ERROR trajCost SmolCrank1 requires field type linear or cubic.\n");
      exit(-1);
    }
    solnGrid = dynamic_cast<const Piecewise1d*>(fieldList[fieldSel.get(2)]);
    if (solnGrid == NULL) {
      fprintf(stderr,"ERROR trajCost SmolCrank1 requires field type linear or cubic.\n");
      exit(-1);
    }

    //std::numeric_limits<double>::infinity();

    T = 1;
    X = 0;

    // The solution will be stored in these objects.
    soln.n = solnGrid->length();
    soln.dx = solnGrid->getDr();
    soln.dxSq = soln.dx*soln.dx;
    soln.x0 = solnGrid->getR0();
    for (int i = 0; i < soln.n; i++) all.add(i);
    
    // The solution will be stored in these objects.
    solver = new CrankNicholsonSolver1(solnGrid, timestep, 1.0/beta);

    // For each node, make an array of solutions.
    // We are assuming that the dt of all events are the same.
    solnProb = new PiecewiseCubic*[soln.n];
    for (int i = 0; i < soln.n; i++) {
      solnProb[i] = new PiecewiseCubic(*solnGrid);
      solnProb[i]->zero();
    }

    // Check the event delT.
    delT = event[0].del[T];
    for (int e = 0; e < eventNum; e++) {
      double d = fabs(delT - event[e].del[T])/delT;
       if (d > 0.1) {
	 fprintf(stderr,"event 0 del[T] %.15g event %d del[T] %.15g\n", delT, e, event[e].del[T]);
	 fprintf(stderr,"ERROR trajCost SmolCrank1 requires all delta time to be the same.\n");
	 exit(-1);
       }
       if (d > 1e-6) {
	 fprintf(stderr,"event 0 del[T] %.15g event %d del[T] %.15g\n", delT, e, event[e].del[T]);
	 fprintf(stderr,"Warning trajCost SmolCrank1 requires all delta time to be the same.\n");
       }
    }
    steps = int(ceil(delT/timestep));

    // // Validate the diffusivity and force.
    // if (diffuse->length() != force->length() || diffuse->length() != solnGrid->length()) {
    //   fprintf(stderr,"diffuse->length() %d force->length() %d solnGrid->length() %d\n", diffuse->length(), force->length(), solnGrid->length());
    //   fprintf(stderr,"ERROR trajCost SmolCrank1 requires all fields to have the same number of nodes\n");
    //   exit(-1);
    // }

    updateLocal();
    cloneLast();
  }

  ~TrajSmolCrank1() {
    for (int i = 0; i < soln.n; i++) delete[] solnProb[i];
    delete[] solnProb;
    delete solver;
  }


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

  // A region around the node defined by the user.
  IndexList nodeRegion(int home) const {
    if (maxHop >= soln.n) return all;
    IndexList ret;

    if (solnGrid->getPeriodic()) {
      ret.add(home);
      for (int i = 1; i < maxHop; i++) {
	ret.add(solnGrid->wrapNode(home-i));
	ret.add(solnGrid->wrapNode(home+i));
      }
    } else {
      int i0 = home-maxHop;
      int i1 = home+maxHop;
      if (i0 < 0) i0 = 0;
      if (i1 >= soln.n) i1 = soln.n-1;
      for (int i = i0; i <= i1; i++) ret.add(i);
    }
    return ret;
  }

  // Get the cost for an event by using the Smoluchowski solution
  // over the diffusivity and force fields.
  double eventCost(int e) {
    double x0 = event[e].var[X];
    int node0 = solnGrid->nearestNode(x0);
    double x1 = x0 + event[e].del[X];

    // We have stored a solution vector for each node.
    double p1 = solnProb[node0]->computeVal(x1);
    // Check for invalid probabilities.
    if (p1 <= 0.0 || p1 != p1) return std::numeric_limits<double>::infinity();

    // The cost
    return -log(p1);
  }

  void solve(int i) {
    solver->init(solnProb[i], solnProb[i]->getNodePos(i), soln.dx);
    solver->solve(solnProb[i], steps, diffuse, force, NULL);
    solver->conserveProb(solnProb[i]);

    int j = solver->positive(solnProb[i]);
    if (j >= 0) {
      fprintf(stderr,"Warning. Negative probablity %g at position %d. Invalidating this solution.\n", solnProb[i]->get(j), j);

      //PiecewiseCubic g(soln.x0, soln.x0 + soln.dx*soln.n, solnProb[i], soln.n, false);
      //g.write("neg.dat");
      //diffuse->write("neg.diffuse");

      // Invalidate the whole solution.
      solnProb[i]->zero();
      solnProb[i]->shift(-1);
    }
  }

private:

  ////////////////////////////////////////////////////////////////////
    public:
  // We overload deltaCost so that we can precompute the Smoluchowski solution.
  virtual double deltaCost(const TrialMove& trialMove) {
    int node = refField->nearestNode(trialMove.pos.x);
    IndexList region = nodeRegion(node);
    
    // Solve the Smoluchowski equation for each node in the region.
#pragma omp parallel for
    for (int i = 0; i < region.length(); i++) {
      int i0 = region.get(i);
      solve(i0);
    }

    return TrajCostComputer::deltaCost(trialMove, region);
  }

  // We overload cost so that we can precompute the Smoluchowski solution.
  virtual double calcCost() {
    // Solve the Smoluchowski equation for all nodes.
#pragma omp parallel for
    for (int i = 0; i < soln.n; i++) solve(i);

    return TrajCostComputer::calcCost();
  }

  // We overload revert so that we can recompute the Smoluchowski solution.
  virtual void revert(const TrialMove& trialMove) {
    int node = refField->nearestNode(trialMove.pos.x);
    IndexList region = nodeRegion(node);

    // Re-solve the Smoluchowski equation for each node in the region.
#pragma omp parallel for
    for (int i = 0; i < region.length(); i++) {
      int i0 = region.get(i);
      solve(i0);
    }

    TrajCostComputer::revert(trialMove, region);
  }


  // Set the local costs.
  virtual double updateLocal() {
    // Solve the Smoluchowski equation for all nodes.
#pragma omp parallel for
    for (int i = 0; i < soln.n; i++) solve(i);

    return TrajCostComputer::updateLocal();
  }
};

#endif
