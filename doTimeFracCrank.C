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

#include "useful.H"
#include "PiecewiseZero.H"
#include "PiecewiseCubic.H"
#include "TimeFracCrankSolver.H"
#include <omp.h>

int main(int argc, char* argv[]) {
  if ( argc != 9 ) {
    printf("Usage: %s inFile diffuseFile forceFile alphaFile timestep time kT outFile\n", argv[0]);
    return 0;
  }
  const char* inFile = argv[1];
  const char* diffuseFile = argv[2];
  const char* forceFile = argv[3];
  const char* alphaFile = argv[4];
  double timestep = strtod(argv[5], NULL);
  double tim = strtod(argv[6], NULL);
  double kT = strtod(argv[7], NULL);
  const char* outFile = argv[argc-1];

  printf("\ndoTimeFracCrank\n");
#pragma omp parallel
  {
    printf("Threads: %d\n", omp_get_num_threads());
  }

  bool periodic = false;

  // Load the diffusivity and force.
  PiecewiseCubic diffuse(diffuseFile, periodic);
  PiecewiseCubic force(forceFile, periodic);
  PiecewiseCubic alpha(alphaFile, periodic);

  // Load the solution domain.
  PiecewiseZero init(inFile, periodic);
  const int n = init.length();
  double* prob = new double[n];
  for (int i = 0; i < n; i++) prob[i] = init.get(i);

  // Prepare to do the finite difference.
  int steps = int(ceil(tim/timestep));
  if (steps == 0) {
    fprintf(stderr, "No steps!\n");
    exit(-1);
  }
  printf("steps %d\n", steps);

  // The main object.
  TimeFracCrankSolver solver(&init, timestep, kT);

  // Initialize the clock.
  //clock_t clockInit = clock();
  double clockInit = omp_get_wtime();

  solver.solve(prob, steps, &diffuse, &force, NULL, &alpha);
  init.reset(prob, init.getPeriodic());
  double scale = solver.conserveProb(init);
  printf("scale %g\n", scale);
  init.write(outFile);
  printf("value1 %g\n", prob[n/2]);

  //printf("\nRun time: %.4g s\n", double(clock()-clockInit)/CLOCKS_PER_SEC);
  double runTime = omp_get_wtime() - clockInit;
  printf("Run time: %.10g s\n", runTime);

  delete[] prob;
  return 0;
}
