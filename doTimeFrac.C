///////////////////////////////////////////////////////////////////////
// Author: Jeff Comer <jeffcomer at gmail>

#include <cstdio>
#include <cmath>

#include "useful.H"
#include "PiecewiseZero.H"
#include "PiecewiseCubic.H"
#include "TimeFracCrankSolver.H"

int main(int argc, char* argv[]) {
  if ( argc != 10 ) {
    printf("Usage: %s inFile diffuseFile forceFile alphaFile timestep time kT outPeriod outPrefix\n", argv[0]);
    return 0;
  }
  const char* inFile = argv[1];
  const char* diffuseFile = argv[2];
  const char* forceFile = argv[3];
  const char* alphaFile = argv[4];
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
  PiecewiseCubic alpha(alphaFile, periodic);

  // Load the solution domain.
  PiecewiseZero init(inFile, periodic);
  const int n = init.length();
  double* prob = new double[n];
  for (int i = 0; i < n; i++) prob[i] = init.get(i);

  // Prepare to do the finite difference.
  int steps = int(ceil(tim/timestep));
  if (outPeriod <= 0) outPeriod = steps;
  int cycles = steps/outPeriod;
  printf("steps %d\n", steps);
  printf("cycles %d\n", cycles);

  // The main object.
  TimeFracCrankSolver solver(&init, timestep, kT);

  char outFile[256];
  for (int c = 0; c < cycles; c++) {
    // Write the output.
    init.reset(prob, init.getPeriodic());
    solver.conserveProb(init);
    snprintf(outFile, 256, "%s.%d.dat", outPrefix, c);
    init.write(outFile);

    solver.solve(prob, outPeriod, &diffuse, &force, NULL, &alpha);
  }

  init.reset(prob, init.getPeriodic());
  solver.conserveProb(init);
  snprintf(outFile, 256, "%s.final.dat", outPrefix);
  init.write(outFile);
  printf("value1 %g\n", prob[n/2]);

  delete[] prob;
  return 0;
}

