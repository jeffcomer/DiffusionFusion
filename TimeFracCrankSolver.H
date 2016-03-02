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
#ifndef TIMEFRACCRANKSOLVER_H
#define TIMEFRACCRANKSOLVER_H

#include "Piecewise1d.H"
#include <omp.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_gamma.h>

class TimeFracCrankSolver {
private:
  int n;
  double dx, x0;
  bool periodic;
  double timestep;
  double beta;
  static const int jM = 20;

public:
  TimeFracCrankSolver(const Piecewise1d* refField, double timestep0, double kT) {
    n = refField->length();
    dx = refField->getDr();
    x0 = refField->getR0();
    periodic = refField->getPeriodic();
    timestep = timestep0;
    beta = 1.0/kT;
  }

  int length() const { return n; }

  // Compute 1/(x Gamma(x)) safely.
  static double x_gamma_x_recip(double x) {
    // https://en.wikipedia.org/wiki/Reciprocal_gamma_function
    if (x < 1e-6) return (1.0 + M_EULER*x);
    return gsl_sf_gammainv(x)/x;
  }

  static double memorySum(const int i, const int s, double al, const gsl_vector** solList) {
    double sum = 0.0;
    for (int j = 2; j <= s; j++) {
      double omega = pow(j,1.0-al) - pow(j-1,1.0-al);
      sum += omega*(gsl_vector_get(solList[s-j+1], i) - gsl_vector_get(solList[s-j], i));
    }
    return sum;
  }

  // Nice attempt, but isn't accurate.
  static double memorySumApprox(const int i, const int s, double al, const gsl_vector** solList, double& longSum) {
    double sum = 0.0;
    int stop = (s>jM)?jM:s;
    for (int j = 2; j <= stop; j++) {
      double omega = pow(j,1.0-al) - pow(j-1,1.0-al);
      sum += omega*(gsl_vector_get(solList[s-j+1], i) - gsl_vector_get(solList[s-j], i));
    }
    
    if (s > jM) {
      double newTerm = (1.0 - al)*pow(jM+1,-al)*(gsl_vector_get(solList[s-jM+1], i) - gsl_vector_get(solList[s-jM-1], i));
      // long sum and update it at the same time.
      longSum = newTerm + (1.0 - 2.0*al/(s + jM))*longSum;
    
      return sum + longSum;
    } else {
      longSum = 0.0;
      return sum;
    }
  }

  inline double sigmaAlpha(double a) const {
    return 2.0*x_gamma_x_recip(1.0-a)*pow(timestep,-a);
  }
		   

  void init(double* prob, int node) const {
    // Zero the initial distribution.
    for (int i = 0; i < n; i++) prob[i] = 0.0;
    // Add the delta function at the initial point.
    prob[node] = 1.0/dx;
  }

  // Based on N. H. SWEILAM, M. M. KHADER, A. M. S. MAHDY. CRANK-NICOLSON FINITE DIFFERENCE METHOD FOR SOLVING TIME-FRACTIONAL DIFFUSION EQUATION (2012) Journal of Fractional Calculus and Applications 2(2):1-9
  // Equation 8.
  // Multiply gamma -> r and sigma_a,k -> 2*dif*sigma_a,k
  void solve(double* prob, int steps, const Piecewise1d* diffuse, const Piecewise1d* force, const Piecewise1d* bias, const Piecewise1d* alpha) const {
    // Coefficients contain the diffusion and drift effects.
    double ca[n];
    double cb[n];
    double r[n];
    double alphaL[n];
    //double longSum[n];
    double sigma[n];

    for (int i = 0; i < n; i++) {
      // Get the diffusivity and force.
      double dif, gradDif;
      double frc, gradFrc;
      double fb, gradFb;
      double x = x0 + i*dx;
      diffuse->computeValGrad(x, dif, gradDif);
      force->computeValGrad(x, frc, gradFrc);
      alphaL[i] = alpha->computeVal(x);
      sigma[i] = sigmaAlpha(alphaL[i]);

      if (bias != NULL) {
       	bias->computeValGrad(x, fb, gradFb);
       	frc += fb;
       	gradFrc += gradFb;
       }

      /////////////////////////////////////////////////////////////////
      // We are solving the differential equation for p(x,t)
      // p^alpha_t = D*p_xx + (D_x - beta*D*F)*p_x - beta*(D*F_x + D_x*F)*p
      // Here, subscripts represent partial derivatives.
      // D(x): the diffusivity
      // F(x,t) = -w_x(x) + fBias(x,t): the total force
      // beta: 1/(kT)
      //
      // We write this equation as
      // p^alpha_t = ca*p + cb*p_x + D*p_xx
      // Here r differs from Numerical Methods Using Matlab 4th Ed. Mathews and Fink, pp. 561-562
      // because we have divided out timestep.
      ca[i] = -beta*(dif*gradFrc + gradDif*frc);
      cb[i] = 0.5*(gradDif - beta*dif*frc)/dx;
      r[i] = dif/(dx*dx);

      // The long timescale sum starts at zero.
      //longSum[i] = 0.0;
    }

    // Boundary conditions.
    // Left
    double frcA = force->computeVal(x0);
    if (bias != NULL) frcA += bias->computeVal(x0);
    // A dimensionless variable related to the drift at the left boundary.
    double driftA = beta*frcA*dx; 
    // Right
    double x1 = x0 + (n-1)*dx;
    double frcB = force->computeVal(x1);
    if (bias != NULL) frcB += bias->computeVal(x1);
    // A dimensionless variable related to the drift at the right boundary.
    double driftB = beta*frcB*dx;
    
    // Allocate linear system vectors.
    // We have a tridiagonal matrix.
    // See: https://www.gnu.org/software/gsl/manual/html_node/Tridiagonal-Systems.html
    gsl_vector* d = gsl_vector_alloc(n); // Matrix diagonal
    gsl_vector* e = gsl_vector_alloc(n-1); // Band to the right of diagonal
    gsl_vector* f = gsl_vector_alloc(n-1); // Band to the left of diagonal
    gsl_vector* b = gsl_vector_alloc(n); // the right side of Mx = b

    gsl_vector* solList[steps+1]; // the solution "x" at each step.
    for (int s = 0; s <= steps; s++) solList[s] = gsl_vector_alloc(n);

    // Fill in the matrix.
    // Robin boundary conditions (zero flux).
    // Assume alpha = 1 at the boundaries.
    // Left: Use the O(dx^2) forward finite difference.
    {
      int i = 0;
      double corrA = (cb[i]-r[i])/(3.0+2.0*driftA);
      gsl_vector_set(e, i, -(cb[i] + r[i] + corrA));
    }
    // Right: Use the O(dx^2) backward finite difference.
    {
      int i = n-1;
      double corrB = -(cb[i]+r[i])/(3.0-2.0*driftB);
      gsl_vector_set(f, i-1, cb[i] - (r[i] + corrB));
    }
    // Interior nodes using a Crank-Nicholson approach.
    for (int i = 1; i < n-1; i++) {
      //gsl_vector_set(d, i, 2.0 + 2.0*r[i] - ca[i]);
      // The diagonal changes with the timestep due to fractional derivative.
      gsl_vector_set(f, i-1, cb[i] - r[i]);
      gsl_vector_set(e, i, -(cb[i] + r[i]));
    }
    
    // Absorbing at left.
    //gsl_vector_set(d, 0, 1.0);
    //gsl_vector_set(e, 0, 0.0);
    //gsl_vector_set(b, 0, 0.0);
    // Absorbing at right.
    //gsl_vector_set(f, n-2, 0.0);
    //gsl_vector_set(d, n-1, 1.0);
    //gsl_vector_set(b, n-1, 0.0);

    // Initialize the solution vector.
    for (int i = 0; i < n; i++) gsl_vector_set(solList[0], i, prob[i]);

    for (int s = 1; s <= steps; s++) {
      // Set the pointers
      gsl_vector* sol0 = solList[s-1]; // Last solution.
      gsl_vector* sol = solList[s]; // This solution.

      // Fill in the vector b.
      // Reflecting at left.
      {
	int i = 0;
	double um = (4.0*gsl_vector_get(sol0,i) - gsl_vector_get(sol0,i+1))/(3.0+2.0*driftA);
	
	//double sum = memorySumApprox(i, s, alphaL[i], (const gsl_vector**)solList, longSum[i]);
	double sum = memorySum(i, s, alphaL[i], (const gsl_vector**)solList);
	double corrA = (cb[i]-r[i])/(3.0+2.0*driftA);

	gsl_vector_set(d, i, sigma[i] + 2.0*r[i] - ca[i] + 4.0*corrA);
	gsl_vector_set(b, i,
		       (sigma[i] - 2.0*r[i] + ca[i])*gsl_vector_get(sol0,i)
		       + (r[i] - cb[i])*um
		       + (r[i] + cb[i])*gsl_vector_get(sol0,i+1)
		       - sigma[i]*sum);
      }
      // Reflecting at right.
      {
	int i = n-1;
	double up = (4.0*gsl_vector_get(sol0,i) - gsl_vector_get(sol0,i-1))/(3.0-2.0*driftB);

	//double sum = memorySumApprox(i, s, alphaL[i], (const gsl_vector**)solList, longSum[i]);
	double sum = memorySum(i, s, alphaL[i], (const gsl_vector**)solList);
	double corrB = -(cb[i]+r[i])/(3.0-2.0*driftB);

	gsl_vector_set(d, i, sigma[i] + 2.0*r[i] - ca[i] + 4.0*corrB);
	gsl_vector_set(b, i,
		       (sigma[i] - 2.0*r[i] + ca[i])*gsl_vector_get(sol0,i)
		       + (r[i] - cb[i])*gsl_vector_get(sol0,i-1)
		       + (r[i] + cb[i])*up
		       - sigma[i]*sum);
      }

      // Interior
#pragma omp parallel for
      for (int i = 1; i < n-1; i++) {
      double sum = memorySum(i, s, alphaL[i], (const gsl_vector**)solList);
      //double sum = memorySumApprox(i, s, alphaL[i], (const gsl_vector**)solList, longSum[i]);

      gsl_vector_set(d, i, sigma[i] + 2.0*r[i] - ca[i]);
      gsl_vector_set(b, i,
	(sigma[i] - 2.0*r[i] + ca[i])*gsl_vector_get(sol0, i)
	+ (r[i] - cb[i])*gsl_vector_get(sol0, i-1)
	+ (r[i] + cb[i])*gsl_vector_get(sol0, i+1)
	- sigma[i]*sum);
      }

      // Solve the system of linear equations.
      gsl_linalg_solve_tridiag(d, e, f, b, sol);
    }
	  
    // Copy the results into prob.
    for (int i = 0; i < n; i++) prob[i] = gsl_vector_get(solList[steps], i);

    // Deallocate everything.
    gsl_vector_free(d);
    gsl_vector_free(e);
    gsl_vector_free(f);
    gsl_vector_free(b);
    for (int s = 0; s <= steps; s++) gsl_vector_free(solList[s]);
  }

  double conserveProb(double* sol) const {
    double sum = 0.0;
    for (int i = 0; i < n; i++) sum += sol[i];
    double scale = 1.0/(sum*dx);
    //printf("scale %g\n", sum*dx);
    for (int i = 0; i < n; i++) sol[i] *= scale;
    return scale;
  }

  double conserveProb(Piecewise1d& sol) const {
    double sum = 0.0;
    for (int i = 0; i < sol.length(); i++) sum += sol.get(i);
    double scale = 1.0/(sum*sol.getDr());
    //printf("scale %g\n", sum*sol.getDr());
    sol.scale(scale);
    return scale;
  }

  int positive(const double* sol) const {
    for (int i = 0; i < n; i++) {
      if (sol[i] < -0.5/dx) return i;
    }
    return -1;
  }

  ~TimeFracCrankSolver() {
  }
};

#endif
