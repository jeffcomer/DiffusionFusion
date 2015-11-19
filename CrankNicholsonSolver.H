///////////////////////////////////////////////////////////////////////
// Interface for a Smoluchowski equation solver
// Author: Jeff Comer <jeffcomer at gmail>
#ifndef CRANKNICHOLSONSOLVER_H
#define CRANKNICHOLSONSOLVER_H

#include "Piecewise1d.H"
#include "SolutionStruct.H"
#include <omp.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>

class CrankNicholsonSolver {
private:
  int n;
  double dx, x0;
  bool periodic;
  double timestep;
  double invTimestep;
  double beta;

public:
  CrankNicholsonSolver(const Piecewise1d* refField, double timestep0, double kT) {
    n = refField->length();
    dx = refField->getDr();
    x0 = refField->getR0();
    periodic = refField->getPeriodic();
    timestep = timestep0;
    invTimestep = 1.0/timestep;
    beta = 1.0/kT;
  }

  void init(double* prob, int node) const {
    // Zero the initial distribution.
    for (int i = 0; i < n; i++) prob[i] = 0.0;
    // Add the delta function at the initial point.
    prob[node] = 1.0/dx;
  }

  void solve(double* prob, int steps, const Piecewise1d* diffuse, const Piecewise1d* force, const Piecewise1d* bias) const {
    // Coefficients contain the diffusion and drift effects.
    double ca[n];
    double cb[n];
    double r[n];
    for (int i = 0; i < n; i++) {
      // Get the diffusivity and force.
      double dif, gradDif;
      double frc, gradFrc;
      double fb, gradFb;
      double x = x0 + i*dx;
      diffuse->computeValGrad(x, dif, gradDif);
      force->computeValGrad(x, frc, gradFrc);
      // dif = diffuse->get(i);
      // gradDif = diffuse->getGrad(i);
      // frc = force->get(i);
      // gradFrc = force->getGrad(i);

      if (bias != NULL) {
       	bias->computeValGrad(x, fb, gradFb);
       	frc += fb;
       	gradFrc += gradFb;
       }

      // Calculate the coefficients.
      ca[i] = -beta*(dif*gradFrc + gradDif*frc)*timestep;
      cb[i] = 0.5*(gradDif - beta*dif*frc)/dx*timestep;
      r[i] = dif*timestep/(dx*dx);
    }

    // Boundary conditions.
    //double fudge = 1.17; // Why?
    // Left
    double frcA = force->computeVal(x0);
    if (bias != NULL) frcA += bias->computeVal(x0);
    double driftA = beta*frcA*dx;
    // Right
    double x1 = x0 + (n-1)*dx;
    double frcB = force->computeVal(x1);
    if (bias != NULL) frcB += bias->computeVal(x1);
    double driftB = beta*frcB*dx;
    
    // Allocate linear system vectors.
    gsl_vector* d = gsl_vector_alloc(n);
    gsl_vector* e = gsl_vector_alloc(n-1);
    gsl_vector* f = gsl_vector_alloc(n-1);
    gsl_vector* b = gsl_vector_alloc(n);
    gsl_vector* sol = gsl_vector_alloc(n);
    gsl_vector* sol0 = gsl_vector_alloc(n);

    // Intialize the solution vector.
    for (int i = 0; i < n; i++) gsl_vector_set(sol, i, prob[i]);

    for (int s = 0; s < steps; s++) {
      // Swap the pointers.
      gsl_vector* tmp = sol0;
      sol0 = sol;
      sol = tmp;
      
      // Robin boundary conditions (zero flux).
      // Use the O(dx^2) forward finite difference.
      // Reflecting at left.
      {
	int i = 0;
	double corr = (cb[i]-r[i])/(3.0+2.0*driftA);
	gsl_vector_set(d, i, 2.0 + 2.0*r[i] - ca[i] + 4.0*corr);
	gsl_vector_set(e, i, -(cb[i] + r[i] + corr));
	double um = (4.0*gsl_vector_get(sol0,i) - gsl_vector_get(sol0,i+1))/(3.0+2.0*driftA);
	gsl_vector_set(b, i,
		       (2.0 - 2.0*r[i] + ca[i])*gsl_vector_get(sol0,i)
		       + (r[i] - cb[i])*um
		       + (r[i] + cb[i])*gsl_vector_get(sol0,i+1));
      }

      // Reflecting at right.
      // Use the O(dx^2) backward finite difference.
      // The finite difference formulas change signs; hence, 3.0-2.0*driftB.
      {
	int i = n-1;
	double corr = -(cb[i]+r[i])/(3.0-2.0*driftB);
	gsl_vector_set(d, i, 2.0 + 2.0*r[i] - ca[i] + 4.0*corr);
	gsl_vector_set(f, i-1, cb[i] - (r[i] + corr));
	double up = (4.0*gsl_vector_get(sol0,i) - gsl_vector_get(sol0,i-1))/(3.0-2.0*driftB);
	gsl_vector_set(b, i,
		       (2.0 - 2.0*r[i] + ca[i])*gsl_vector_get(sol0,i)
		       + (r[i] - cb[i])*gsl_vector_get(sol0,i-1)
		       + (r[i] + cb[i])*up);
      }

      // Absorbing at left.
      //gsl_vector_set(d, 0, 1.0);
      //gsl_vector_set(e, 0, 0.0);
      //gsl_vector_set(b, 0, 0.0);

      // Absorbing at right.
      //gsl_vector_set(f, n-2, 0.0);
      //gsl_vector_set(d, n-1, 1.0);
      //gsl_vector_set(b, n-1, 0.0);

      // Interior using a Crank-Nicholson approach.
      // See: Numerical Methods Using Matlab 4th Ed. Mathews and Fink, pp. 561-562
      // See: https://www.gnu.org/software/gsl/manual/html_node/Tridiagonal-Systems.html
      for (int i = 1; i < n-1; i++) {
	gsl_vector_set(d, i, 2.0 + 2.0*r[i] - ca[i]);
	gsl_vector_set(f, i-1, cb[i] - r[i]);
	gsl_vector_set(e, i, -(cb[i] + r[i]));
	gsl_vector_set(b, i,
		       (2.0 - 2.0*r[i] + ca[i])*gsl_vector_get(sol0, i)
		       + (r[i] - cb[i])*gsl_vector_get(sol0, i-1)
		       + (r[i] + cb[i])*gsl_vector_get(sol0, i+1));
      }

      // Solve the system of linear equations.
      gsl_linalg_solve_tridiag(d, e, f, b, sol);
    }

    // Conserve the probability.
    conserveProb0(sol);
	  
    // Copy the results into prob.
    for (int i = 0; i < n; i++) prob[i] = gsl_vector_get(sol, i);

    // Deallocate everything.
    gsl_vector_free(d);
    gsl_vector_free(e);
    gsl_vector_free(f);
    gsl_vector_free(b);
    gsl_vector_free(sol);
    gsl_vector_free(sol0);
  }

  void conserveProb0(gsl_vector* sol) const {
    double sum = 0.0;
    for (int i = 0; i < n; i++) sum += gsl_vector_get(sol, i);
    double scale = 1.0/(sum*dx);
    printf("scale %g\n", sum*dx);
    for (int i = 0; i < n; i++)
      gsl_vector_set(sol, i, scale*gsl_vector_get(sol, i));
  }

  ~CrankNicholsonSolver() {
  }
};

#endif
