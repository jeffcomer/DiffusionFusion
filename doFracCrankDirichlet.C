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
// Author: Jeff Comer <jeffcomer at gmail>

#include <cstdio>
#include <cmath>
#include <ctime>

#include "useful.H"
#include "PiecewiseZero.H"
#include "PiecewiseCubic.H"
#include "TimeFracCrankSolverDirichlet.H"

void calcCurrent(PiecewiseCubic& curr, const PiecewiseCubic& prob, const PiecewiseCubic& diffuse, const PiecewiseCubic& force, double kT) {
  for (int i = 0; i < curr.length(); i++) {
    double x = curr.nodePos(i);
    double conc, gradConc;

    prob.computeValGrad(x,conc,gradConc);
    double dif = diffuse.computeVal(x);
    double frc = force.computeVal(x);
    
    curr.set(i, -dif*gradConc + frc*dif*conc/kT);
  }
}

void calcResistivity(PiecewiseCubic& resisty, const PiecewiseCubic& prob, const PiecewiseCubic& diffuse, const PiecewiseCubic& force, double kT) {
  for (int i = 0; i < resisty.length(); i++) {
    double x = resisty.nodePos(i);
    double conc, gradConc;

    prob.computeValGrad(x,conc,gradConc);
    double dif = diffuse.computeVal(x);
    double frc = force.computeVal(x);
    
    double curr = -dif*gradConc + frc*dif*conc/kT;
    resisty.set(i, -gradConc/curr);
  }
}


int main(int argc, char* argv[]) {
  if ( argc != 10 ) {
    printf("Usage: %s inFile diffuseFile forceFile alphaFile timestep time kT perturbation outPrefix\n", argv[0]);
    return 0;
  }
  const char* inFile = argv[1];
  const char* diffuseFile = argv[2];
  const char* forceFile = argv[3];
  const char* alphaFile = argv[4];
  double timestep = strtod(argv[5], NULL);
  double tim = strtod(argv[6], NULL);
  double kT = strtod(argv[7], NULL);
  double perturbation = strtod(argv[8], NULL);
  const char* outPrefix = argv[argc-1];

  printf("perturbation %g\n", perturbation);
  bool periodic = false;

  // Load the diffusivity and force.
  PiecewiseCubic diffuse(diffuseFile, periodic);
  PiecewiseCubic force(forceFile, periodic);
  PiecewiseCubic alpha(alphaFile, periodic);

  // Load the solution domain.
  PiecewiseCubic init(inFile, periodic);
  const int n = init.length();
  double* prob = new double[n];
  for (int i = 0; i < n; i++) prob[i] = init.get(i);
  // For storing the current.
  PiecewiseCubic curr(init);

  // Prepare to do the finite difference.
  int steps = int(ceil(tim/timestep));
  if (steps == 0) {
    fprintf(stderr, "No steps!\n");
    exit(-1);
  }
  printf("steps %d\n", steps);

  // The main object.
  TimeFracCrankSolverDirichlet solver(&init, timestep, kT, 1.0, 1.0-perturbation);

  // Initialize the clock.
  clock_t clockInit = clock();
  solver.solve(prob, steps, &diffuse, &force, NULL, &alpha);

  char outFile[256];
  // Write the output.
  init.reset(prob, init.getPeriodic());
  snprintf(outFile, 256, "%s.prob", outPrefix);
  init.write(outFile);

  snprintf(outFile, 256, "%s.curr", outPrefix);
  calcCurrent(curr, init, diffuse, force, kT);
  curr.write(outFile);

  snprintf(outFile, 256, "%s.resisty", outPrefix);
  calcResistivity(curr, init, diffuse, force, kT);
  curr.write(outFile);
  printf("value1 %g\n", prob[n/2]);
  
  printf("\nRun time: %.4g s\n", double(clock()-clockInit)/CLOCKS_PER_SEC);

  delete[] prob;
  return 0;
}

