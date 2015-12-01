///////////////////////////////////////////////////////////////////////
// Author: Jeff Comer <jeffcomer at gmail>

#include <cstdio>
#include <cmath>

#include "useful.H"
#include "PiecewiseZero.H"
#include "PiecewiseCubic.H"
#include "CrankNicholsonSolver.H"
#include <omp.h>

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

  #pragma omp parallel
  {
    printf("Threads: %d\n", omp_get_num_threads());
  }

  // Prepare to do the finite difference.
  int steps = int(ceil(tim/timestep));
  int cycles = steps/outPeriod;
  printf("steps %d\n", steps);
  printf("cycles %d\n", cycles);

  // The main object.
  CrankNicholsonSolver solver(&init, timestep, kT);

  // Initialize the clock.
  double clockInit = omp_get_wtime();
  
  char outFile[256];
  for (int c = 0; c < cycles; c++) {
    // Write the output.
    init.reset(prob, init.getPeriodic());
    solver.conserveProb(init);
    snprintf(outFile, 256, "%s.%d.dat", outPrefix, c);
    init.write(outFile);

    solver.solve(prob, outPeriod, &diffuse, &force, &bias);
  }

  snprintf(outFile, 256, "%s.final.dat", outPrefix);
  init.write(outFile);
  
  double runTime = omp_get_wtime() - clockInit;
  printf("\nRun time: %.10g s\n", runTime);

  delete[] prob;
  return 0;
}

