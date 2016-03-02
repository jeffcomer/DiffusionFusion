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
///////////////////////////////////////////////////////////////////////
// Interface for a Smoluchowski equation solver
// Author: Jeff Comer <jeffcomer at gmail>
#ifndef BENESISMOLSOLVER_H
#define BENESISMOLSOLVER_H

#include "PiecewiseCubic.H"
#include "SmoluchowskiSolver.H"

class BenesiSmolSolver : public SmoluchowskiSolver {
private:
  double beta; // 1/(kT)
  int steps; // the number of subintervals of the propagation time
  const int binNum;
  // We swap pointers during the iteration.
  PiecewiseCubic* pbCurr;
  PiecewiseCubic* pbNext;
  double* binCen;
  double* binX;
  double binDx;
  double binHalfDx;

public:
  BenesiSmolSolver(const Field& binDim) : binNum(binDim.length()) {

    const Piecewise1d* pb = dynamic_cast<const Piecewise1d*>(&binDim);
    if (pb == NULL) {
       fprintf(stderr,"ERROR BenesiSmolSolver::BenesiSmolSolver requires a bin of type Piecewise1d.\n");
       exit(-1);
    }
    pbCurr = new PiecewiseCubic(*pb);
    pbNext = new PiecewiseCubic(*pb);

    binDx = pbCurr->getDr();
    binHalfDx = 0.5*binDx;
    binCen = new double[binNum];
    binX = new double[binNum];
    for (int i = 0; i < binNum; i++) {
      binCen[i] = pbCurr->getNodeCenter(i);
      binX[i] = pbCurr->getNodePos(i);
    }
  }

  virtual void setParameters(double beta0, int steps0) {
    beta = fabs(beta0);
    steps = abs(steps0);
  }

  virtual ~BenesiSmolSolver() {
    delete pbCurr;
    delete pbNext;
    delete[] binCen;
    delete[] binX;
  }

  int wrapNode(int j) const {
    if (pbCurr->getPeriodic()) return pbCurr->wrapNode(j);
    // Reflecting boundary conditions.
    if (j < 0) return 0;
    else if (j >= binNum) return binNum-1;
    return j;
  }

  bool writeCurrent(const char* fileName) const {
    return pbCurr->write(fileName);
  }

  virtual double solveProb(double xi, double xf, double tau, const Piecewise1d* diffuse, const Piecewise1d* forceSys, const Piecewise1d* forceBias, double& hopMax) {
    // The timestep.
    double dt = tau/steps;

    // Get the start and end bins.
    int indI = wrapNode(pbCurr->nearestNode(xi));
    int indF = wrapNode(pbCurr->nearestNode(xf));

    // Start with a "delta function" at the bin nearest xi.
    // These are probability densities.
    pbCurr->zero();
    pbCurr->set(indI, 1.0/binDx);
    pbNext->zero();
    int regionNum = 0;

    // Run the finite difference steps.
    for (int s = 0; s < steps; s++) {
      // Clear the relevent part of the bin for the next step.
      pbNext->zero();
      regionNum = 0;

      /// Write a file.
      //char outFile[256];
      //snprintf(outFile, 256, "output/bin.%d.dat", s);
      //pbCurr->write(outFile);

      // Run over all nodes.
      for (int i = 0; i < binNum; i++) {
	// Ignore untouched nodes.
	if (pbCurr->get(i) == 0.0) continue;
	regionNum++;

	// Get the diffusivity.
	double dif, gradDif;
	diffuse->computeValGrad(binCen[i],dif,gradDif);
	double frc = forceSys->computeVal(binCen[i]);
	// We add the bias force if there is one.
	if (forceBias != NULL) frc += forceSys->computeZeroOrder(binCen[i]);
	  
	// The probability in the bin has drifted and expanded by:
	double drift = (beta*frc*dif + gradDif)*dt;
	double sigma = sqrt(2.0*dif*dt);
	double halfWidth = binHalfDx + sigma;
	// The probability amplitude:
	double probAmpNext = pbCurr->get(i)*binHalfDx/halfWidth;

	// These nodes bracket the transformed prob.
	double x0 = binCen[i] + drift - halfWidth;
	double x1 = binCen[i] + drift + halfWidth;
	int i0 = pbCurr->getNode(x0);
	int i1 = pbCurr->getNode(x1);
	// Wrap the nodes, but keep the old logical versions, because
	// they allow us to calculate the remaining pieces.
	int iw0 = wrapNode(i0);
	int iw1 = wrapNode(i1);

	if (iw0 == iw1) {
	  // All of the probability goes into one node.
	  // This shouldn't happen very often--likely at the boundaries.
	  pbNext->set(i0, pbNext->get(iw0) + pbCurr->get(i)*binDx);
	} else {
	  // The contribution to the first and last nodes.
	  // getNodePos will return logical positions for unwrapped node
	  // indices, which is what we want.
	  double p0 = probAmpNext*(pbCurr->getNodePos(i0+1) - x0)/binDx;
	  double p1 = probAmpNext*(x1 - pbCurr->getNodePos(i1))/binDx;
	  pbNext->set(i0, pbNext->get(i0) + p0);
	  pbNext->set(i1, pbNext->get(i1) + p1);

	  // Nodes in between.
	  if (!pbCurr->getPeriodic()) {
	    // We just add the amplitude over the whole node.
	    // j is guaranteed to be valid running from iw0 to iw1.
	    for (int j = iw0+1; j<=iw1-1; j++) {
	      pbNext->set(j, pbNext->get(j) + probAmpNext);
	    }
	  } else {
	    // We just add the amplitude over the whole node.
	    // We run over the logical nodes and wrap as necessary.
	    for (int j = i0+1; j<=i1-1; j++) {
	      int jw = pbCurr->wrapNode(j);
	      pbNext->set(jw, pbNext->get(jw) + probAmpNext);
	    } 
	  } // end in-between node (next) periodicity conditional
	  //printf("s %d i %d i0 %d i1 %d prob0 %g prob %g\n", s, i, i0, i1, pbCurr->get(i), pSum);
	} // conditional for multiple nodes (next)
      } // end loop over relevant nodes (current)
      
      // Cycle the pointers.
      PiecewiseCubic* tmp = pbCurr;
      pbCurr = pbNext;
      pbNext = tmp;
      //printf("regionNum %d\n", regionNum);
      double probDen = 0.0;
      for (int i = 0; i < binNum; i++) probDen += pbCurr->get(i);
      //printf("prob %g\n", probDen*binDx);
    }

    if (0.5*binDx*regionNum > hopMax) hopMax = 0.5*binDx*regionNum;
    // Return the probability in the desired bin.
    return pbCurr->get(indF)*binDx;
  }

  
};

#endif
