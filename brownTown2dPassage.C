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

#include <ctime>
#include "useful.H"
#include "PiecewiseBicubic.H"
#include "RandomGsl.H"

double reflectPolar(double x, double eps) {
  // Handle very long jumps.
  if (x < -180.0 + eps) x = -180.0 + eps;
  else if (x > 180.0 - eps) x = 180.0 -eps;
  
  // Reflecting boundary.
  if (x < 0) x = -x;
  else if (x > 180.0) x = 360.0 - x;
  return x;
}

int main(int argc, char* argv[]) {
  if (argc != 13) {
    printf("Usage: %s pmfFile diffuseZFile diffuseThFile z0 goalZ th0 dt kT steps outPeriod seed0 outFile\n", argv[0]);
    exit(0);
  }
  const char* outFile = argv[argc-1];
  PiecewiseBicubic pmf(argv[1], false);  
  PiecewiseBicubic diffuseZ(argv[2], false);  
  PiecewiseBicubic diffuseTh(argv[3], false);
  double z0 = strtod(argv[4],NULL);
  double goalZ = strtod(argv[5],NULL);
  double th0 = strtod(argv[6],NULL);
  double dt = strtod(argv[7],NULL);
  double kT = strtod(argv[8],NULL);
  long int steps = atol(argv[9]);
  long int seed0 = atol(argv[10]);
  int outPeriod = atoi(argv[11]);

  printf("PMF `%s' %d nodes\n", argv[1], pmf.length());
  printf("diffuseZ `%s' %d nodes\n", argv[2], diffuseZ.length());
  printf("diffuseTh `%s' %d nodes\n", argv[3], diffuseTh.length());
  printf("z0 %g goalZ %g th0 %g dt %g kT %g steps %ld outPeriod %d\n", z0, goalZ, th0, dt, kT, steps, outPeriod);

  double beta = 1.0/kT;
  long seed = (unsigned int)time((time_t *)NULL) + seed0*seed0*seed0;
  Random rando(seed);

  // Open the file.
  FILE* out = fopen(outFile,"w");
  if (out == NULL) {
    printf("ERROR Couldn't open file %s.\n",outFile);
    return false;
  }
  
  fprintf(out, "%.15g %.15g %.15g\n", 0.0, z0, th0);
  double lastZ = z0;
  double lastTh = th0;
  double currZ = lastZ;
  double currTh = lastTh;
  long int s;
  for (s = 1; s <= steps; s++) {
    double t = dt*s;

    // Adding in the effective force would be double counting, since
    // it's in the PMF we get from colvars
    
    // Calculate force and diffusivity.
    double forZ = -pmf.computeGradX(lastZ,lastTh);
    double forTh = -pmf.computeGradY(lastZ,lastTh);
    double difZ = diffuseZ.computeVal(lastZ,lastTh);
    double gradDifZ = diffuseZ.computeGradX(lastZ,lastTh);
    double difTh = diffuseTh.computeVal(lastZ,lastTh);
    double gradDifTh = diffuseTh.computeGradY(lastZ,lastTh);

    currZ = lastZ + (beta*forZ*difZ + gradDifZ)*dt + sqrt(2.0*difZ*dt)*rando.gaussian();
    double sigTh = sqrt(2.0*difTh*dt);
    currTh = lastTh + (beta*forTh*difZ + gradDifTh)*dt + sigTh*rando.gaussian();

    // Boundary conditions.
    currTh = reflectPolar(currTh, sigTh);
    currZ = pmf.wrapX(currZ);

    if (s % outPeriod == 0) 
      fprintf(out, "%.15g %.15g %.15g\n", t, currZ, currTh);
    // Prepare for the next iteration.
    lastZ = currZ;
    lastTh = currTh;

    // Stop if we have reached the goal.
    if (goalZ > z0 && currZ > goalZ) break;
    else if (goalZ < z0 && currZ < goalZ) break;
  }

  if (s <= steps) printf ("Reached goal at t=%g\n", dt*s);
  else printf("Did not reach goal.\n");


  fclose(out);

  return 0;
}
