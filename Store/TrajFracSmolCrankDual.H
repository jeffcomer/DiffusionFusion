// Author: Jeff Comer <jeffcomer at gmail>
#ifndef TRAJFRACSMOLCRANKDUAL_H
#define TRAJFRACSMOLCRANKDUAL_H

#include "useful.H"
#include "Field.H"
#include "PiecewiseLinear.H"
#include "PiecewiseCubic.H"
#include "PiecewiseZero.H"
#include "Event.H"
#include "TrajCostComputer.H"
#include "SolutionStruct.H"
#include "TimeFracCrankSolver.H"

class TrajFracSmolCrankDual : public TrajCostComputer {
private:
  const Piecewise1d* diffuse;
  const Piecewise1d* force;
  const Piecewise1d* alpha;
  int T, X;
  const Piecewise1d* refField;
  double timestep, maxHop;
  
  SolutionStruct soln;
  double** solnProb;
  double** solnProb1;
  double delT, delT1;
  int steps, steps1;
  IndexList all;
  TimeFracCrankSolver* solver;
  
public:
  TrajFracSmolCrankDual(const Field** ppf, IndexList fieldSel0, const Event* event0, int eventNum0, double kt, int leastLocal, double timestep0, double maxHop0) :
    TrajCostComputer(ppf, fieldSel0, event0, eventNum0, kt, leastLocal, 2), timestep(timestep0), maxHop(maxHop0) {

    if (tcFieldNum != 3) {
      fprintf(stderr, "ERROR trajCost FracSmolCrankDual takes three fields: (0) diffusivity (1) force (2) alpha\n");
      exit(-1);
    }

    diffuse = dynamic_cast<const Piecewise1d*>(fieldList[fieldSel.get(0)]);
    if (diffuse == NULL) {
      fprintf(stderr,"ERROR trajCost FracSmolCrankDual requires field type linear or cubic.\n");
      exit(-1);
    }
    force = dynamic_cast<const Piecewise1d*>(fieldList[fieldSel.get(1)]);
    if (force == NULL) {
      fprintf(stderr,"ERROR trajCost FracSmolCrankDual requires field type linear or cubic.\n");
      exit(-1);
    }
    alpha = dynamic_cast<const Piecewise1d*>(fieldList[fieldSel.get(2)]);
    if (alpha == NULL) {
      fprintf(stderr,"ERROR trajCost FracSmolCrankDual requires field type linear or cubic.\n");
      exit(-1);
    }

    T = 1;
    X = 0;

    refField = dynamic_cast<const Piecewise1d*>(leastLocalField); // defined in TrajCostComputer
    if (refField == NULL) {
      fprintf(stderr,"ERROR trajCost FracSmolCrankDual requires leastLocalField type linear or cubic.\n");
      exit(-1);
    }

    // The solution will be stored in these objects.
    soln.n = refField->length();
    soln.dx = refField->getDr();
    soln.dxSq = soln.dx*soln.dx;
    soln.x0 = refField->getR0();
    for (int i = 0; i < soln.n; i++) all.add(i);
    
    // The solution will be stored in these objects.
    solver = new TimeFracCrankSolver(refField, timestep, 1.0/beta);

    // For each node, make an array of solutions.
    // We assume that we have only two distinct dt values.
    solnProb = new double*[soln.n];
    solnProb1 = new double*[soln.n];
    for (int i = 0; i < soln.n; i++) {
      solnProb[i] = new double[soln.n];
      solnProb1[i] = new double[soln.n];
      for (int j = 0; j < soln.n; j++) {
	solnProb[i][j] = 0.0;
	solnProb1[i][j] = 0.0;
      }
    }

    // Check the event delT.
    delT = event[0].del[T];
    delT1 = -1;
    for (int e = 0; e < eventNum; e++) {
      double d1;
      double d = fabs(delT - event[e].del[T])/delT;
      if (d > 1e-6 && delT1 < 0.0) {
	delT1 = event[e].del[T];
      } else {
	d1 = fabs(delT1 - event[e].del[T])/delT;
	if (d > 1e-6 && d1 > 1e-6) {
	  fprintf(stderr,"delT %.15g delT1 %.15g event %d del[T] %.15g\n", delT, delT1, e, event[e].del[T]);
	  fprintf(stderr,"ERROR trajCost FracSmolCrankDual requires no more than two distinct delta time values.\n");
	  exit(-1);
	}
      }
    }
    fprintf(stdout, "delT %.15g delT1 %.15g\n", delT, delT1);

    steps = int(ceil(delT/timestep));
    steps1 = int(ceil(delT1/timestep));
    updateLocal();
    cloneLast();
  }

  ~TrajFracSmolCrankDual() {
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
    case 2:
      return String("alpha");
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

    double cost;
    if (fabs(delT - event[e].del[T]) < fabs(delT1 - event[e].del[T])) {
      // This event has delT.
      if (solnProb[node0][node1] <= 0.0 || solnProb[node0][node1] != solnProb[node0][node1]) {
	return std::numeric_limits<double>::infinity();
      }
      cost = -log(solnProb[node0][node1]);
    } else {
      // This event has delT1.
      if (solnProb1[node0][node1] <= 0.0 || solnProb1[node0][node1] != solnProb1[node0][node1]) {
	return std::numeric_limits<double>::infinity();
      }
      cost = -log(solnProb1[node0][node1]);
    }

    return cost;
  }

  void solve(int i) {
    int j;

    solver->init(solnProb[i], i);
    solver->solve(solnProb[i], steps, diffuse, force, NULL, alpha);
    solver->conserveProb(solnProb[i]);

    j = solver->positive(solnProb[i]);
    if (j >= 0) {
      fprintf(stdout,"Warning. Negative probablity %g at position %d. Invalidating this solution.\n", solnProb[i][j], j);
      // Invalidate the whole solution.
      for (int k = 0; k < soln.n; k++) solnProb[i][k] = -1.0;
    }


    solver->init(solnProb1[i], i);
    solver->solve(solnProb1[i], steps, diffuse, force, NULL, alpha);
    solver->conserveProb(solnProb1[i]);

    j = solver->positive(solnProb1[i]);
    if (j >= 0) {
      fprintf(stdout,"Warning. Negative probablity %g at position %d. Invalidating this solution.\n", solnProb1[i][j], j);
      // Invalidate the whole solution.
      for (int k = 0; k < soln.n; k++) solnProb1[i][k] = -1.0;
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
