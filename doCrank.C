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
#include "CrankNicholsonSolver.H"

int main(int argc, char* argv[]) {
  if ( argc != 10 ) {
    //printf("Usage: %s inFile diffuseFile forceFile biasFile timestep time kT periodic outPeriod outPrefix\n", argv[0]);
    printf("Usage: %s inFile diffuseFile forceFile biasFile timestep time kT outPeriod outPrefix\n", argv[0]);
    return 0;
  }
  const char* inFile = argv[1];
  const char* diffuseFile = argv[2];
  const char* forceFile = argv[3];
  const char* biasFile = argv[4];
  double timestep = strtod(argv[5], NULL);
  double tim = strtod(argv[6], NULL);
  double kT = strtod(argv[7], NULL);
  //int periodic = atoi(argv[8]);
  int outPeriod = atoi(argv[argc-2]);
  const char* outPrefix = argv[argc-1];

  bool periodic = false;

  // Load the diffusivity and force.
  PiecewiseCubic diffuse(diffuseFile, periodic);
  PiecewiseCubic force(forceFile, periodic);
  PiecewiseCubic bias(biasFile, periodic);

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
  if (outPeriod <= 0) outPeriod = steps;
  int cycles = steps/outPeriod;
  printf("steps %d\n", steps);
  printf("cycles %d\n", cycles);

  // The main object.
  CrankNicholsonSolver solver(&init, timestep, kT);

  // Initialize the clock.
  clock_t clockInit = clock();

  char outFile[256];
  for (int c = 0; c < cycles; c++) {
    // Write the output.
    init.reset(prob, init.getPeriodic());
    solver.conserveProb(init);
    snprintf(outFile, 256, "%s.%d.dat", outPrefix, c);
    init.write(outFile);

    solver.solve(prob, outPeriod, &diffuse, &force, &bias);
  }

  init.reset(prob, init.getPeriodic());
  solver.conserveProb(init);
  snprintf(outFile, 256, "%s.final.dat", outPrefix);
  init.write(outFile);
  printf("value1 %g\n", prob[n/2]);
  
  printf("\nRun time: %.4g s\n", double(clock()-clockInit)/CLOCKS_PER_SEC);

  delete[] prob;
  return 0;
}

