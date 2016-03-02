// Author: Jeff Comer <jeffcomer at gmail>
#ifndef TRAJFDSMOLUCHOWSKI_H
#define TRAJFDSMOLUCHOWSKI_H

#include "useful.H"
#include "Field.H"
#include "PiecewiseLinear.H"
#include "PiecewiseCubic.H"
#include "PiecewiseZero.H"
#include "Event.H"
#include "TrajCostComputer.H"
#include "SolutionStruct.H"

class TrajFdSmoluchowski : public TrajCostComputer {
private:
  const Piecewise1d* diffuse;
  const Piecewise1d* force;
  const Piecewise1d** biasFieldList;
  int T, X, D;
  // Maximum node hop made by events in each node.
  int* nodeHopMax;
  const Piecewise1d* refField;
  double timestep;

  SolutionStruct soln;
  double* solnProb;
  double* arDif;
  double* arGradDif;
  double* arFrc;
  double* arGradFrc;
  
public:
  TrajFdSmoluchowski(const Field** ppf, IndexList fieldSel0, const Event* event0, int eventNum0, double kt, int leastLocal, const Piecewise1d** biasFieldList0, double timestep0) :
    TrajCostComputer(ppf, fieldSel0, event0, eventNum0, kt, leastLocal, 3), biasFieldList(biasFieldList0), nodeHopMax(NULL), timestep(timestep0) {

    if (tcFieldNum != 2) {
      fprintf(stderr, "ERROR trajCost smoluchowski takes two fields: (0) diffusivity (1) force\n");
      exit(-1);
    }

    diffuse = dynamic_cast<const Piecewise1d*>(fieldList[fieldSel.get(0)]);
    if (diffuse == NULL) {
      fprintf(stderr,"ERROR trajCost smoluchowski requires field type linear or cubic.\n");
      exit(-1);
    }
    force = dynamic_cast<const Piecewise1d*>(fieldList[fieldSel.get(1)]);
    if (force == NULL) {
      fprintf(stderr,"ERROR trajCost smoluchowski requires field type linear or cubic.\n");
      exit(-1);
    }

    //std::numeric_limits<double>::infinity();

    T = 1;
    X = 0;
    D = 0;

    refField = dynamic_cast<const Piecewise1d*>(leastLocalField); // defined in TrajCostComputer
    if (refField == NULL) {
      fprintf(stderr,"ERROR trajCost smoluchowski requires leastLocalField type linear or cubic.\n");
      exit(-1);
    }

    // The solution will be stored in this objects.
    soln.n = refField->length();
    soln.dx = refField->getDr();
    soln.dxSq = soln.dx*soln.dx;
    soln.x0 = refField->getR0();
    solnProb = new double[soln.n];
    arDif = new double[soln.n];
    arGradDif = new double[soln.n];
    arFrc = new double[soln.n];
    arGradFrc = new double[soln.n];
    // We can interpolate the diffusivity and force.
    for (int i = 0; i < soln.n; i++) {
      double x = refField->getNodePos(i);
      diffuse->computeValGrad(x,arDif[i],arGradDif[i]);
      force->computeValGrad(x,arFrc[i],arGradFrc[i]);
    }

    // For each node, we look for the neighboring nodes affected by events beginning in this node.
    nodeHopMax = new int[leastLocalNodes];
    for (int n = 0; n < leastLocalNodes; n++) nodeHopMax[n] = 3;
    
    updateLocal();
    updateHopMax();
    cloneLast();
  }

  ~TrajFdSmoluchowski() {
    delete[] solnProb;
    delete[] arDif;
    delete[] arGradDif;
    delete[] arFrc;
    delete[] arGradFrc;
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

  IndexList nodeRegion(int home) const {
    IndexList region;
    int n0 = home - nodeHopMax[home];
    int n1 = home + nodeHopMax[home];

    if (leastLocalField->getPeriodic()) {
      // Account for wrapping and make sure nodes are unique.
      for (int n = n0; n <= n1; n++) {
	int nw = refField->wrapNode(n);
	if (nw == n || nw < n0) region.add(n);
      }
    } else {
      if (n0 < 0) n0 = 0;
      if (n1 >= leastLocalNodes) n1 = leastLocalNodes-1;

      for (int n = n0; n <= n1; n++) region.add(n);
    }
    return region;
  }

  // Estimate the maximum hop distance for BD from each node.
  void updateHopMax() {
    for (int n = 0; n < leastLocalNodes; n++) {
      int ne = local[n].events.length();

      // Find the biggest hop from each node.
      int hopMax = 0;
      for (int le = 0; le < ne; le++) {
	int e = local[n].events.get(le);
	// Since we do this without omp parallel for, the distribution should be
	// in solnProb.
	eventCost(e);
	
	double x0 = event[e].var[X];
	int homeNode = refField->nearestNode(x0);
	double sumProb = soln.dx*solnProb[homeNode];
	int i0 = homeNode;
	int i1 = homeNode;
	int count = 0;
	// Find how far we have to go to contain 0.95 prob.
	if (refField->getPeriodic()) {
	  while (sumProb < 0.95) {
	    i0 = refField->wrapNode(i0-1);
	    i1 = refField->wrapNode(i1+1);
	    // We quit if these reach the same point.
	    if (i0 == i1) break;

	    sumProb += soln.dx*solnProb[i0];
	    sumProb += soln.dx*solnProb[i1];
	    count++;
	  }
	} else {
	  while (sumProb < 0.95) {
	    i0--;
	    i1++;
	    // We quit if we've reached both ends.
	    if (i0 < 0 && i1 >= refField->length()) break;

	    if (i0 >= 0 && i0 < refField->length())
	      sumProb += soln.dx*solnProb[i0];
	    if (i1 <= 0 && i1 < refField->length())
	      sumProb += soln.dx*solnProb[i1];
	    count++;
	  }
	}

	// Count is a good measure of the maximum radius.
	if (count > hopMax) hopMax = count;
      }
      
      // The +3 is because the change at i affects the
      // interpolation at i+2.
      int oldHopMax = nodeHopMax[n];
      nodeHopMax[n] = hopMax + 3;
      //#ifdef FUSION_DEBUG
	if (abs(oldHopMax - nodeHopMax[n]) > 2) printf("nodeHopMax node %d changed from %d to %d\n", n, oldHopMax, nodeHopMax[n]);
      //#endif
    }
  }

  // Get the cost for an event by running the Smoluchowski over the diffusivity and force fields.
  double eventCost(int e) {
    double x0 = event[e].var[X];
    double x1 = x0 + event[e].del[X];

    double* arBias = NULL;
    double* arGradBias = NULL;
    if (event[e].bias >= 0) {
      const Piecewise1d* biasField = biasFieldList[event[e].bias];
      arBias = new double[soln.n];
      arGradBias = new double[soln.n];
      for (int i = 0; i < soln.n; i++) {
	double x = refField->getNodePos(i);
	biasField->computeValGrad(x, arBias[i], arGradBias[i]);
      }
    }

    // Each call to event cost has its own solution vector, so that
    // they can be run in parallel.
    // We have thread one store the solution vector permanently, so
    // that other functions can use the solution.
    double* prob;
    if (omp_get_thread_num() == 0) prob = solnProb;
    else prob = new double[soln.n];
    
    // Zero the initial distribution.
    for (int i = 0; i < soln.n; i++) prob[i] = 0.0;
    // Add the delta function at the initial point.
    int node0 = refField->nearestNode(x0);
    int node1 = refField->nearestNode(x1);
    prob[node0] = 1.0/soln.dx;

    // Solve the Smoluchowski equation on these distributions.
    if (event[e].bias < 0) smoluchowskiSolve(prob, event[e].del[T]);
    else smoluchowskiSolveBias(prob, event[e].del[T], arBias, arGradBias);
    if (prob[node1] == 0.0 || prob[node1] == -0.0) {
      fprintf(stderr,"node0 %d x0 %g node %d x1 %g prob %g\n", node0, x0, node1, x1, prob[node1]);
      fprintf(stderr,"nodes %d x0 %g dx %g\n", soln.n, soln.x0, soln.dx);
      fprintf(stderr,"ERROR trajCost smoluchowski: Probability of arrival is zero. Check -timestep.\n");
      exit(-1);
      return 0.0;
    }
    double cost = -log(prob[node1]);

    // double dist[2];
    // dist[0] = refField->getR0();
    // dist[1] = refField->getR1();
    // PiecewiseCubic test(dist, prob, soln.n, refField->getPeriodic());
    // test.write("test_prob.dat");
    //printf("eventCost x0 %g x1 %g prob %g cost %g\n", x0, x1, prob[node1], cost);
    // exit(0);
    
    if (prob != solnProb) delete[] prob;
    if (arBias != NULL) {
      delete[] arBias;
      delete[] arGradBias;
    }
    return cost;
  }

private:
  // Iterate the Smoluchowski equation on a give probability distribution.
  // The diffusivity and force arrays must be set.
  void smoluchowskiSolve(double* prob, double dt) {
    double* prob0 = new double[soln.n];
    double* p0 = prob0;
    double* p = prob;

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

    //printf("i %d dif %g frc %g\n", i, arDif[i], arFrc[i]);
    return arDif[i]*ddpi + (arGradDif[i] - arDif[i]*arFrc[i])*dpi
      - (arDif[i]*arGradFrc[i] + arGradDif[i]*arFrc[i])*pc;
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
  // Do everything again allowing for a bias. Not very elegant.

  // Iterate the Smoluchowski equation on a give probability distribution.
  // The diffusivity and force arrays must be set.
  void smoluchowskiSolveBias(double* prob, double dt, const double* arBias, const double* arGradBias) {
    double* prob0 = new double[soln.n];
    double* p0 = prob0;
    double* p = prob;
    
    // Steps is always even so the solution is always in prob.
    int steps = 2*int(0.5*dt/timestep);
    if (steps < 4) {
      fprintf(stderr,"ERROR trajCost smoluchowski: Smoluchowski solver steps < 4. Reduce -timestep.\n");
      exit(-1);
    }

    // Iterate the solution.
    if (refField->getPeriodic()) {
      for (int s = 0; s < steps; s++) {
	// Swap the pointers.
	double* tmp = p0;
	p0 = p;
	p = tmp;

	// Do the finite difference calculation.
	finiteDiffPeriodicBias(p0, p, arBias, arGradBias);
      }
    } else {
      for (int s = 0; s < steps; s++) {
	// Swap the pointers.
	double* tmp = p0;
	p0 = p;
	p = tmp;

	// Do the finite difference calculation.
	finiteDiffBias(p0, p, arBias, arGradBias);
      }
    }

    delete[] prob0;
  }

  // The left-hand side of the differential equation.
  double slopeBias(int i, double pm1, double pc, double pp1, const double* arBias, const double* arGradBias) {
    // 3-point numerical derivatives.
    double dpi = 0.5*(pp1-pm1)/soln.dx;
    double ddpi = (pm1 - 2.0*pc + pp1)/soln.dxSq;

    double frc = arFrc[i] + arBias[i];
    double gradFrc = arGradFrc[i] + arGradBias[i];

    return arDif[i]*ddpi + (arGradDif[i] - arDif[i]*frc)*dpi
      - (arDif[i]*gradFrc + arGradDif[i]*frc)*pc;
  }

  void finiteDiffPeriodicBias(const double* p0, double* p, const double* arBias, const double* arGradBias) {
    // Absorbing boundary conditions at left.
    {
      int i = 0;
      double pm1 = p0[soln.n-1];
      double pi = p0[i];
      double pp1 = p0[i+1];

      p[i] = pi + slopeBias(i,pm1,pi,pp1,arBias,arGradBias)*timestep;
    }

    // Absorbing boundary conditions at right.
    {
      int i = soln.n-1;
      double pm1 = p0[i-1];
      double pi = p0[i];
      double pp1 = p0[0];

      p[i] = pi + slopeBias(i,pm1,pi,pp1,arBias,arGradBias)*timestep;
    }

    // The interior.
    for (int i = 1; i < soln.n-1; i++) {
      double pm1 = p0[i-1];
      double pi = p0[i];
      double pp1 = p0[i+1];

      p[i] = pi + slopeBias(i,pm1,pi,pp1,arBias,arGradBias)*timestep;
    }
  }

  void finiteDiffBias(const double* p0, double* p, const double* arBias, const double* arGradBias) {
    // Absorbing boundary conditions at left.
    {
      int i = 0;
      double pm1 = 0.0;
      double pi = p0[i];
      double pp1 = p0[i+1];

      p[i] = pi + slopeBias(i,pm1,pi,pp1,arBias,arGradBias)*timestep;
    }

    // Absorbing boundary conditions at right.
    {
      int i = soln.n-1;
      double pm1 = p0[i-1];
      double pi = p0[i];
      double pp1 = 0.0;

      p[i] = pi + slopeBias(i,pm1,pi,pp1,arBias,arGradBias)*timestep;
    }

    // The interior.
    for (int i = 1; i < soln.n-1; i++) {
      double pm1 = p0[i-1];
      double pi = p0[i];
      double pp1 = p0[i+1];

      p[i] = pi + slopeBias(i,pm1,pi,pp1,arBias,arGradBias)*timestep;
    }
  }

  ////////////////////////////////////////////////////////////////////
  // Calculate the cost over all nodes.
    public:
  virtual double calcCost() {
    // We can interpolate the diffusivity and force for this cycle.
    for (int i = 0; i < soln.n; i++) {
      double x = refField->getNodePos(i);
      diffuse->computeValGrad(x,arDif[i],arGradDif[i]);
      force->computeValGrad(x,arFrc[i],arGradFrc[i]);
    }

    long double cost = 0.0;
#pragma omp parallel for reduction(+:cost)
    for (int e = 0; e < eventNum; e++) cost += eventCost(e);
    return cost;
  }

  // We overload deltaCost, revert because we can no longer rely on just nearby neighbors.
  virtual double deltaCost(const TrialMove& trialMove) {
    printf("deltaCost\n");
    // We can interpolate the diffusivity and force for this cycle.
    for (int i = 0; i < soln.n; i++) {
      double x = refField->getNodePos(i);
      diffuse->computeValGrad(x,arDif[i],arGradDif[i]);
      force->computeValGrad(x,arFrc[i],arGradFrc[i]);
    }

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
