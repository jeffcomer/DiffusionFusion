///////////////////////////////////////////////////////////////////////
// Interface for a Smoluchowski equation solver
// Author: Jeff Comer <jeffcomer at gmail>
#ifndef SMOLUCHOWSKISOLVER_H
#define SMOLUCHOWSKISOLVER_H

#include "Piecewise1d.H"
#include "SolutionStruct.H"
#include <omp.h>

struct DiffusiveParameters {
public:
  DiffusiveParameters(const Piecewise1d* d, const Piecewise1d* f, const Piecewise1d* b):
    diffuse(d), force(f), bias(b) {
  }

  const Piecewise1d* diffuse;
  const Piecewise1d* force;
  const Piecewise1d* bias;
};

class SmoluchowskiSolver {
private:
  int n;
  double dx, dxSq, x0;
  bool periodic;
  double timestep;
  double beta;

public:
  SmoluchowskiSolver(const Piecewise1d* refField, double timestep0, double kT) {
    n = refField->length();
    dx = refField->getDr();
    dxSq = dx*dx;
    x0 = refField->getR0();
    periodic = refField->getPeriodic();
    timestep = timestep0;
    beta = 1.0/kT;
  }

  void init(double* prob, int node) {
    // Zero the initial distribution.
    for (int i = 0; i < n; i++) prob[i] = 0.0;
    // Add the delta function at the initial point.
    prob[node] = 1.0/dx;
  }

  void solve(double* prob, int steps, const Piecewise1d* diffuse, const Piecewise1d* force, const Piecewise1d* bias) {
    double prob0[n];
    double* p0 = prob0;
    double* p = prob;

    DiffusiveParameters parm(diffuse,force,bias);

    // Iterate the solution.
    if (periodic) {
      for (int s = 0; s < steps; s++) {
	// Swap the pointers.
	double* tmp = p0;
	p0 = p;
	p = tmp;

	// Do the finite difference calculation.
	finiteDiffPeriodic(p0, p, parm);
      }
    } else {
      for (int s = 0; s < steps; s++) {
	// Swap the pointers.
	double* tmp = p0;
	p0 = p;
	p = tmp;

	// Do the finite difference calculation.
	finiteDiff3(p0, p, parm);
      }
    }

    // If we ran an odd number steps, copy the contents of p into prob.
    if (p != prob) {
      for (int i = 0; i < n; i++) prob[i] = p[i];
    }
  }

  void finiteDiffPeriodic(const double* p0, double* p, DiffusiveParameters& parm) {
    // Periodic boundary conditions at left.
    {
      int i = 0;
      double pm1 = p0[n-1];
      double pi = p0[i];
      double pp1 = p0[i+1];

      p[i] = pi + slope(i,pm1,pi,pp1,parm)*timestep;
    }

    // Periodic boundary conditions at right.
    {
      int i = n-1;
      double pm1 = p0[i-1];
      double pi = p0[i];
      double pp1 = p0[0];

      p[i] = pi + slope(i,pm1,pi,pp1,parm)*timestep;
    }

    // The interior.
#pragma omp parallel for
    for (int i = 1; i < n-1; i++) {
      double pm1 = p0[i-1];
      double pi = p0[i];
      double pp1 = p0[i+1];

      p[i] = pi + slope(i,pm1,pi,pp1,parm)*timestep;
    }
  }

  void finiteDiff(const double* p0, double* p, DiffusiveParameters& parm) {
    // Reflecting boundary conditions at left.
    {
      int i = 0;
      double pm1 = p0[i];
      double pi = p0[i];
      double pp1 = p0[i+1];

      p[i] = pi + slope(i,pm1,pi,pp1,parm)*timestep;
    }

    // Reflecting boundary conditions at right.
    {
      int i = n-1;
      double pm1 = p0[i-1];
      double pi = p0[i];
      double pp1 = p0[i];

      p[i] = pi + slope(i,pm1,pi,pp1,parm)*timestep;
    }

    // The interior.
#pragma omp parallel for
    for (int i = 1; i < n-1; i++) {
      double pm1 = p0[i-1];
      double pi = p0[i];
      double pp1 = p0[i+1];

      p[i] = pi + slope(i,pm1,pi,pp1,parm)*timestep;
    }
  }

  void finiteDiff1(const double* p0, double* p, DiffusiveParameters& parm) {
    // The interior.
#pragma omp parallel for
    for (int i = 1; i < n-1; i++) {
      double pm1 = p0[i-1];
      double pi = p0[i];
      double pp1 = p0[i+1];

      p[i] = pi + slope(i,pm1,pi,pp1,parm)*timestep;
    }

    // Reflecting boundary conditions.
    // Found from no-flux condition.
    p[0] = p[1];
    p[n-1] = p[n-2];
  }

  void finiteDiff2(const double* p0, double* p, DiffusiveParameters& parm) {
    // Reflecting boundary conditions at left.
    {
      int i = 0;
      double pm1 = p0[i+1];
      double pi = p0[i];
      double pp1 = p0[i+1];

      p[i] = pi + slope(i,pm1,pi,pp1,parm)*timestep;
    }

    // Reflecting boundary conditions at right.
    {
      int i = n-1;
      double pm1 = p0[i-1];
      double pi = p0[i];
      double pp1 = p0[i-1];

      p[i] = pi + slope(i,pm1,pi,pp1,parm)*timestep;
    }

    // The interior.
#pragma omp parallel for
    for (int i = 1; i < n-1; i++) {
      double pm1 = p0[i-1];
      double pi = p0[i];
      double pp1 = p0[i+1];

      p[i] = pi + slope(i,pm1,pi,pp1,parm)*timestep;
    }
  }

  void finiteDiffBound(const double* p0, double* p, DiffusiveParameters& parm) {
    // Reflecting boundary conditions at left.
    {
      int i = 0;
      double pi = p0[i];
      p[i] = pi + slopeBound(i,pi,parm)*timestep;
    }

    // Reflecting boundary conditions at right.
    {
      int i = n-1;
      double pi = p0[i];
      p[i] = pi + slopeBound(i,pi,parm)*timestep;
    }

    // The interior.
    for (int i = 1; i < n-1; i++) {
      double pm1 = p0[i-1];
      double pi = p0[i];
      double pp1 = p0[i+1];

      p[i] = pi + slope(i,pm1,pi,pp1,parm)*timestep;
    }
  }

  void finiteDiff3(const double* p0, double* p, DiffusiveParameters& parm) {
    // The interior.
    double sum = 0.0;
#pragma omp parallel for reduction(+:sum)
    for (int i = 1; i < n-1; i++) {
      double pm1 = p0[i-1];
      double pi = p0[i];
      double pp1 = p0[i+1];

      p[i] = pi + slope(i,pm1,pi,pp1,parm)*timestep;
      sum += p[i];
    }

    // Reflecting boundary conditions.
    // Found from no-flux condition.
    sum *= dx;
    //printf("%g\n", sum);
    if (p[1] != 0.0 || p[n-2] != 0.0) {
      double weight = p[1]/(p[1] + p[n-2]);      
      p[0] = (1.0-sum)/dx*weight;
      p[n-1] = (1.0-sum)/dx*(1.0-weight);
    }
  }

  // The left-hand side of the differential equation.
  double slope(int i, double pm1, double pc, double pp1, DiffusiveParameters& parm) {
    // 3-point numerical derivatives.
    double dpi = 0.5*(pp1-pm1)/dx;
    double ddpi = (pm1 - 2.0*pc + pp1)/dxSq;

    double dif, gradDif;
    double frc, gradFrc;
    double bias, gradBias;

    dif = parm.diffuse->get(i);
    gradDif = parm.diffuse->getGrad(i);
    frc = parm.force->get(i);
    gradFrc = parm.force->getGrad(i);

    //parm.diffuse->computeValGrad(x, dif, gradDif);
    //parm.force->computeValGrad(x, frc, gradFrc);

    //double x;
    if (parm.bias != NULL) {
    //x = x0 + i*dx;
    //parm.bias->computeValGrad(x, bias, gradBias);
       bias = parm.bias->get(i);
       gradBias = parm.bias->getGrad(i);
       frc += bias;
       gradFrc += gradBias;
    }

    //return dif*ddpi + (gradDif - dif*frc)*dpi - (dif*gradFrc + gradDif*frc)*pc;
    return dif*ddpi + (gradDif - beta*dif*frc)*dpi - beta*(dif*gradFrc + gradDif*frc)*pc;
  }

  // Reflecting no flux boundary.
  double slopeBound(int i, double pc, DiffusiveParameters& parm) {
    double dif, gradDif;
    double frc, gradFrc;
    double bias, gradBias;

    double x = x0 + i*dx;
    parm.diffuse->computeValGrad(x, dif, gradDif);
    parm.force->computeValGrad(x, frc, gradFrc);

    if (parm.bias != NULL) {
      parm.bias->computeValGrad(x, bias, gradBias);
      frc += bias;
      gradFrc += gradBias;
    }

    return -beta*(dif*gradFrc + gradDif*frc)*pc;
  }

  ~SmoluchowskiSolver() {
  }
};

#endif
