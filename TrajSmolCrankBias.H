// Author: Jeff Comer <jeffcomer at gmail>
#ifndef TRAJSMOLCRANKBIAS_H
#define TRAJSMOLCRANKBIAS_H

#include "useful.H"
#include "Field.H"
#include "PiecewiseLinear.H"
#include "PiecewiseCubic.H"
#include "PiecewiseZero.H"
#include "Event.H"
#include "TrajCostComputer.H"
#include "CrankNicholsonSolver.H"

class TrajSmolCrankBias : public TrajCostComputer {
private:
  const Piecewise1d* diffuse;
  const Piecewise1d* force;
  const Piecewise1d** biasFieldList;
  int T, X;
  const Piecewise1d* refField;
  CrankNicholsonSolver* solver;
  double timestep, maxHop;
  IndexList all;
  
public:
  TrajSmolCrankBias(const Field** ppf, IndexList fieldSel0, const Event* event0, int eventNum0, double kt, int leastLocal, const Piecewise1d** biasFieldList0, double timestep0, double maxHop0) :
    TrajCostComputer(ppf, fieldSel0, event0, eventNum0, kt, leastLocal, 3), biasFieldList(biasFieldList0), timestep(timestep0), maxHop(maxHop0) {

    if (tcFieldNum != 2) {
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
      fprintf(stderr,"ERROR trajCost smolCrank requires leastLocalField type linear or cubic.\n");
      exit(-1);
    }

    // The solution calculator.
    solver = new CrankNicholsonSolver(refField, timestep, 1.0/beta);
    for (int i = 0; i < refField->length(); i++) all.add(i);

    updateLocal();
    cloneLast();
  }

  ~TrajSmolCrankBias() {
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

  // Overload neighbors() to give 
  // a region around the node defined by the user.
  virtual IndexList neighbors(int home) const {
    IndexList ret;

    if (maxHop >= refField->length()) return all;

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
      if (i1 >= refField->length()) i1 = refField->length()-1;
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

    if (node0 < 0 || node0 >= refField->length()) return 0.0;
    if (node1 < 0 || node1 >= refField->length()) return 0.0;

    // Do the solution from scratch for each event. Expensive!
    double solnProb[solver->length()];
    
    // Set the initial condition.
    solver->init(solnProb, node0);
    int steps = int(event[e].del[T]/timestep);
    if (steps < 4) {
      fprintf(stderr,"ERROR trajCost smolCrank: Smoluchowski solver steps < 4. Reduce -timestep.\n");
      exit(-1);
    }

    // Get the bias force field.
    const Piecewise1d* bias;
    if (event[e].bias >= 0) bias = biasFieldList[event[e].bias];
    else bias = NULL;

    // Solve.
    solver->solve(solnProb, steps, diffuse, force, bias);
    solver->conserveProb(solnProb);
    //double scale = solver->conserveProb(solnProb, refField->length());
    /*if (scale < 1e-6 || scale > 10.0) {
      fprintf(stderr,"WARNING trajCost smolCrank: Probability poorly conserved. scale=%g. Increase -biasSmooth or make -timestep or grid size finer.\n", scale);
      return -log(0.0);
      }*/

    if (solnProb[node1] == 0.0 || solnProb[node1] == -0.0 || solnProb[node1] != solnProb[node1]) {
      fprintf(stderr,"node0 %d x0 %g node %d x1 %g prob %g\n", node0, x0, node1, x1, solnProb[node1]);
      fprintf(stderr,"ERROR trajCost smolCrank: Probability of arrival is zero or NaN. Check -timestep.\n");
      return -log(0.0);
    } 
    // We have stored a solution vector for each node.
    double cost = -log(solnProb[node1]);
    //nanCheck(cost,'e',e,event[e].del[X],solnProb[node1],solnProb[node1],solnProb[node1]);
    return cost;
  }
};

#endif
