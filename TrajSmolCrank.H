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
#ifndef TRAJSMOLCRANK_H
#define TRAJSMOLCRANK_H

#include "useful.H"
#include "Field.H"
#include "PiecewiseLinear.H"
#include "PiecewiseCubic.H"
#include "PiecewiseZero.H"
#include "Event.H"
#include "TrajCostComputer.H"
#include "SolutionStruct.H"
#include "CrankNicholsonSolver.H"

class TrajSmolCrank : public TrajCostComputer {
private:
  const Piecewise1d* diffuse;
  const Piecewise1d* force;
  int T, X;
  const Piecewise1d* refField;
  double timestep, maxHop;
  
  SolutionStruct soln;
  double** solnProb;
  double delT;
  int steps;
  IndexList all;
  CrankNicholsonSolver* solver;
  
public:
  TrajSmolCrank(TrajCostDesc& tcd) :
    TrajCostComputer(tcd,2), timestep(tcd.timestep), maxHop(tcd.maxHop) {

    if (fieldSel.length() != 2) {
      fprintf(stderr, "ERROR trajCost smolCrank takes two fields: (0) diffusivity (1) force\n");
      exit(-1);
    }

    diffuse = dynamic_cast<const Piecewise1d*>(fieldList[fieldSel.get(0)]);
    if (diffuse == NULL) {
      fprintf(stderr,"ERROR trajCost smolCrank requires field type linear or cubic.\n");
      exit(-1);
    }
    force = dynamic_cast<const Piecewise1d*>(fieldList[fieldSel.get(1)]);
    if (force == NULL) {
      fprintf(stderr,"ERROR trajCost smolCrank requires field type linear or cubic.\n");
      exit(-1);
    }

    //std::numeric_limits<double>::infinity();

    T = 1;
    X = 0;

    refField = dynamic_cast<const Piecewise1d*>(leastLocalField); // defined in TrajCostComputer
    if (refField == NULL) {
      fprintf(stderr,"ERROR trajCost SmolCrank requires leastLocalField type linear or cubic.\n");
      exit(-1);
    }

    // The solution will be stored in these objects.
    soln.n = refField->length();
    soln.dx = refField->getDr();
    soln.dxSq = soln.dx*soln.dx;
    soln.x0 = refField->getR0();
    for (int i = 0; i < soln.n; i++) all.add(i);
    
    // The solution will be stored in these objects.
    solver = new CrankNicholsonSolver(refField, timestep, 1.0/beta);

    // For each node, make an array of solutions.
    // We are assuming that the dt of all events are the same.
    solnProb = new double*[soln.n];
    for (int i = 0; i < soln.n; i++) {
      solnProb[i] = new double[soln.n];
      for (int j = 0; j < soln.n; j++) solnProb[i][j] = 0.0;
    }

    // Check the event delT.
    delT = event[0].del[T];
    for (int e = eventStart; e <= eventEnd; e++) {
      double d = fabs(delT - event[e].del[T])/delT;
       if (d > 0.1) {
	 fprintf(stderr,"event 0 del[T] %.15g event %d del[T] %.15g\n", delT, e, event[e].del[T]);
	 fprintf(stderr,"ERROR trajCost SmolCrank requires all delta time to be the same.\n");
	 exit(-1);
       }
       if (d > 1e-6) {
	 fprintf(stderr,"event 0 del[T] %.15g event %d del[T] %.15g\n", delT, e, event[e].del[T]);
	 fprintf(stderr,"Warning trajCost SmolCrank requires all delta time to be the same.\n");
       }
    }
    steps = int(ceil(delT/timestep));

    // // Validate the diffusivity and force.
    // if (diffuse->length() != force->length() || diffuse->length() != refField->length()) {
    //   fprintf(stderr,"diffuse->length() %d force->length() %d refField->length() %d\n", diffuse->length(), force->length(), refField->length());
    //   fprintf(stderr,"ERROR trajCost SmolCrank requires all fields to have the same number of nodes\n");
    //   exit(-1);
    // }

    updateLocal();
    cloneLast();
  }

  ~TrajSmolCrank() {
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

    if (refField->getPeriodic()) {
      ret.add(home);
      for (int i = 1; i < maxHop; i++) {
	ret.add(refField->wrapNode(home-i));
	ret.add(refField->wrapNode(home+i));
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
    double x1 = x0 + event[e].del[X];
    int node0 = refField->nearestNode(x0);
    int node1 = refField->nearestNode(x1);

    if (node0 < 0 || node0 >= soln.n) return 0.0;
    if (node1 < 0 || node1 >= soln.n) return 0.0;

    // We have stored a solution vector for each node.
    if (solnProb[node0][node1] == 0.0 || solnProb[node0][node1] == -0.0 || solnProb[node0][node1] != solnProb[node0][node1]) {
      //fprintf(stdout,"node0 %d x0 %g node %d x1 %g prob %g\n", node0, x0, node1, x1, solnProb[node0][node1]);
      //fprintf(stdout,"nodes %d x0 %g dx %g\n", soln.n, soln.x0, soln.dx);
      //fprintf(stderr,"ERROR trajCost smoluchowski: Probability of arrival is zero. Check -timestep.\n");
      //fprintf(stdout,"ERROR trajCost smoluchowski: Probability of arrival is zero. Check -timestep.\n");

      // Write the prob. distribution.
      //PiecewiseCubic out(soln.x0, soln.x0+soln.n*soln.dx, solnProb[node0], soln.n, false);
      //out.write("prob.dat");
      //exit(-1);
      return std::numeric_limits<double>::infinity();
    }
    // We have stored a solution vector for each node.
    double cost = -log(solnProb[node0][node1]);
    return cost;
  }

  void solve(int i) {
    solver->init(solnProb[i], i);
    solver->solve(solnProb[i], steps, diffuse, force, NULL);
    solver->conserveProb(solnProb[i]);

    int j = solver->positive(solnProb[i]);
    if (j >= 0) {
      fprintf(stdout,"Warning. Negative probablity %g at position %d. Invalidating this solution.\n", solnProb[i][j], j);

      //PiecewiseCubic g(soln.x0, soln.x0 + soln.dx*soln.n, solnProb[i], soln.n, false);
      //g.write("neg.dat");
      //diffuse->write("neg.diffuse");

      // Invalidate the whole solution.
      for (int k = 0; k < soln.n; k++) solnProb[i][k] = 0.0;
    }
  }

private:

  ////////////////////////////////////////////////////////////////////
    public:
  // We overload deltaCost so that we can precompute the Smoluchowski solution.
  virtual double deltaCost(const TrialMove& trialMove) {
    IndexList region = nodeRegion(trialMove.node);
    
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
    IndexList region = nodeRegion(trialMove.node);
    
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
