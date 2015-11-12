// Author: Jeff Comer <jeffcomer at gmail>
#ifndef TRAJSIMPLESMOL_H
#define TRAJSIMPLESMOL_H

#include "useful.H"
#include "Field.H"
#include "PiecewiseLinear.H"
#include "PiecewiseCubic.H"
#include "PiecewiseZero.H"
#include "Event.H"
#include "TrajCostComputer.H"
#include "SolutionStruct.H"

class TrajSimpleSmol : public TrajCostComputer {
private:
  const Piecewise1d* diffuse;
  const Piecewise1d* force;
  const Piecewise1d** biasFieldList;
  int T, X, D;
  const Piecewise1d* refField;
  double timestep, maxHop;
  
  SolutionStruct soln;
  double** solnProb;
  double delT;
  IndexList all;
  
public:
  TrajSimpleSmol(const Field** ppf, IndexList fieldSel0, const Event* event0, int eventNum0, double kt, int leastLocal, const Piecewise1d** biasFieldList0, double timestep0, double maxHop0) :
    TrajCostComputer(ppf, fieldSel0, event0, eventNum0, kt, leastLocal, 3), biasFieldList(biasFieldList0), timestep(timestep0), maxHop(maxHop0) {

    if (tcFieldNum != 2) {
      fprintf(stderr, "ERROR trajCost simpleSmol takes two fields: (0) diffusivity (1) force\n");
      exit(-1);
    }

    diffuse = dynamic_cast<const Piecewise1d*>(fieldList[fieldSel.get(0)]);
    if (diffuse == NULL) {
      fprintf(stderr,"ERROR trajCost simpleSmol requires field type linear or cubic.\n");
      exit(-1);
    }
    force = dynamic_cast<const Piecewise1d*>(fieldList[fieldSel.get(1)]);
    if (force == NULL) {
      fprintf(stderr,"ERROR trajCost simpleSmol requires field type linear or cubic.\n");
      exit(-1);
    }

    //std::numeric_limits<double>::infinity();

    T = 1;
    X = 0;
    D = 0;

    refField = dynamic_cast<const Piecewise1d*>(leastLocalField); // defined in TrajCostComputer
    if (refField == NULL) {
      fprintf(stderr,"ERROR trajCost simpleSmol requires leastLocalField type linear or cubic.\n");
      exit(-1);
    }

    // The solution will be stored in these objects.
    soln.n = refField->length();
    soln.dx = refField->getDr();
    soln.dxSq = soln.dx*soln.dx;
    soln.x0 = refField->getR0();
    for (int i = 0; i < soln.n; i++) all.add(i);

    // For each node, make an array of solutions.
    // We are assuming that the dt of all events are the same.
    solnProb = new double*[soln.n];
    for (int i = 0; i < soln.n; i++) {
      solnProb[i] = new double[soln.n];
      for (int j = 0; j < soln.n; j++) solnProb[i][j] = 0.0;
    }

    // Check the event delT.
    delT = event[0].del[T];
    for (int e = 0; e < eventNum; e++) {
      double d = fabs(delT - event[e].del[T])/delT;
       if (d > 0.1) {
	 fprintf(stderr,"event 0 del[T] %.15g event %d del[T] %.15g\n", delT, e, event[e].del[T]);
	 fprintf(stderr,"ERROR trajCost simpleSmol requires all delta time to be the same.\n");
	 exit(-1);
       }
       if (d > 1e-6) {
	 fprintf(stderr,"event 0 del[T] %.15g event %d del[T] %.15g\n", delT, e, event[e].del[T]);
	 fprintf(stderr,"Warning trajCost simpleSmol requires all delta time to be the same.\n");
       }
    }

    // Validate the diffusivity and force.
    if (diffuse->length() != force->length() || diffuse->length() != refField->length()) {
      fprintf(stderr,"diffuse->length() %d force->length() %d refField->length() %d\n", diffuse->length(), force->length(), refField->length());
      fprintf(stderr,"ERROR trajCost simpleSmol requires all fields to have the same number of nodes\n");
      exit(-1);
    }

    updateLocal();
    cloneLast();
  }

  ~TrajSimpleSmol() {
    for (int i = 0; i < soln.n; i++) delete[] solnProb[i];
    delete[] solnProb;
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
    case 2:
      return String("displacement");  
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

    // We have stored a solution vector for each node.
    if (solnProb[node0][node1] == 0.0 || solnProb[node0][node1] == -0.0) {
      fprintf(stderr,"node0 %d x0 %g node %d x1 %g prob %g\n", node0, x0, node1, x1, solnProb[node0][node1]);
      fprintf(stderr,"nodes %d x0 %g dx %g\n", soln.n, soln.x0, soln.dx);
      fprintf(stderr,"ERROR trajCost smoluchowski: Probability of arrival is zero. Check -timestep.\n");
      exit(-1);
      return 0.0;
    }
    // We have stored a solution vector for each node.
    double cost = -log(solnProb[node0][node1]);
    return cost;
  }

private:
  // Iterate the Smoluchowski equation on a give probability distribution.
  // The diffusivity and force arrays must be set.
  void smoluchowskiSolve(int node, double* prob, double dt) {
    double* prob0 = new double[soln.n];
    double* p0 = prob0;
    double* p = prob;

    // Zero the initial distribution.
    for (int i = 0; i < soln.n; i++) p[i] = 0.0;
    // Add the delta function at the initial point.
    p[node] = 1.0/soln.dx;

    // Steps is always even so the solution is always in prob.
    int steps = 2*int(0.5*dt/timestep);
    if (steps < 4) {
      fprintf(stderr,"ERROR trajCost smoluchowski: Smoluchowski solver steps < 4. Reduce -timestep.\n");
      exit(-1);
    }
    //printf("steps %d\n", steps);

    // Iterate the solution.
    if (refField->getPeriodic()) {
     for (int s = 0; s < steps; s++) {
	// Swap the pointers.
	double* tmp = p0;
	p0 = p;
	p = tmp;

	// Do the finite difference calculation.
	finiteDiffPeriodic(p0, p);
      }
     } else {
      for (int s = 0; s < steps; s++) {
	// Swap the pointers.
	double* tmp = p0;
	p0 = p;
	p = tmp;

	// Do the finite difference calculation.
	finiteDiff(p0, p);
      }
    }

    //PiecewiseZero out(p, soln.n, soln.x0, soln.dx, false);
    //out.write("soln.dat");

    delete[] prob0;
  }

  // The left-hand side of the differential equation.
  double slope(int i, double pm1, double pc, double pp1) {
    // 3-point numerical derivatives.
    double dpi = 0.5*(pp1-pm1)/soln.dx;
    double ddpi = (pm1 - 2.0*pc + pp1)/soln.dxSq;

    double dif = diffuse->get(i);
    double gradDif = diffuse->getGrad(i);
    double frc = force->get(i);
    double gradFrc = force->getGrad(i);

    return dif*ddpi + (gradDif - dif*frc)*dpi - (dif*gradFrc + gradDif*frc)*pc;
  }

  void finiteDiffPeriodic(const double* p0, double* p) {
    // Absorbing boundary conditions at left.
    {
      int i = 0;
      double pm1 = p0[soln.n-1];
      double pi = p0[i];
      double pp1 = p0[i+1];

      p[i] = pi + slope(i,pm1,pi,pp1)*timestep;
    }

    // Absorbing boundary conditions at right.
    {
      int i = soln.n-1;
      double pm1 = p0[i-1];
      double pi = p0[i];
      double pp1 = p0[0];

      p[i] = pi + slope(i,pm1,pi,pp1)*timestep;
    }

    // The interior.
    for (int i = 1; i < soln.n-1; i++) {
      double pm1 = p0[i-1];
      double pi = p0[i];
      double pp1 = p0[i+1];

      p[i] = pi + slope(i,pm1,pi,pp1)*timestep;
    }
  }


void finiteDiff(const double* p0, double* p) {
    // Absorbing boundary conditions at left.
    {
      int i = 0;
      double pm1 = 0.0;
      double pi = p0[i];
      double pp1 = p0[i+1];

      p[i] = pi + slope(i,pm1,pi,pp1)*timestep;
    }

    // Absorbing boundary conditions at right.
    {
      int i = soln.n-1;
      double pm1 = p0[i-1];
      double pi = p0[i];
      double pp1 = 0.0;

      p[i] = pi + slope(i,pm1,pi,pp1)*timestep;
    }

    // The interior.
    for (int i = 1; i < soln.n-1; i++) {
      double pm1 = p0[i-1];
      double pi = p0[i];
      double pp1 = p0[i+1];

      p[i] = pi + slope(i,pm1,pi,pp1)*timestep;
    }
  }


  ////////////////////////////////////////////////////////////////////
    public:
  // We overload deltaCost so that we can precompute the Smoluchowski solution.
  virtual double deltaCost(const TrialMove& trialMove) {
    IndexList region = nodeRegion(trialMove.node);
    
    // Solve the Smoluchowski equation for each node in the region.
#pragma omp parallel for
    for (int i = 0; i < region.length(); i++) {
      int i0 = region.get(i);
      smoluchowskiSolve(i0, solnProb[i0], delT);
    }

    return TrajCostComputer::deltaCost(trialMove, region);
  }

  // We overload cost so that we can precompute the Smoluchowski solution.
  virtual double calcCost() {
    // Solve the Smoluchowski equation for all nodes.
#pragma omp parallel for
    for (int i = 0; i < soln.n; i++) smoluchowskiSolve(i, solnProb[i], delT);

    return TrajCostComputer::calcCost();
  }

  // We overload revert so that we can recompute the Smoluchowski solution.
  virtual void revert(const TrialMove& trialMove) {
    IndexList region = nodeRegion(trialMove.node);
    
    // Re-solve the Smoluchowski equation for each node in the region.
#pragma omp parallel for
    for (int i = 0; i < region.length(); i++) {
      int i0 = region.get(i);
      smoluchowskiSolve(i0, solnProb[i0], delT);
    }

    TrajCostComputer::revert(trialMove, region);
  }


  // Set the local costs.
  virtual double updateLocal() {
    // Solve the Smoluchowski equation for all nodes.
#pragma omp parallel for
    for (int i = 0; i < soln.n; i++) smoluchowskiSolve(i, solnProb[i], delT);

    return TrajCostComputer::updateLocal();
  }
};

#endif
