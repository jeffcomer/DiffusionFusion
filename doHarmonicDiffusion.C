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
//
// Based on Schulten, "Nonequilibrium Statistical Mechanics" lecture notes,
// eq. 4.120, 4.121, 4.122

#include <cstdio>
#include <cmath>

#include "useful.H"
#include "PiecewiseZero.H"
#include "PiecewiseCubic.H"

int main(int argc, char* argv[]) {
  if ( argc != 8 ) {
    printf("Usage: %s inFile diffuse spring time beta x0 outFile\n", argv[0]);
    return 0;
  }
  const char* inFile = argv[1];
  double diffuse = strtod(argv[2], NULL);
  double spring = strtod(argv[3], NULL);
  double tim = strtod(argv[4], NULL);
  double beta = strtod(argv[5], NULL);
  double x0 = strtod(argv[6], NULL);
  const char* outFile = argv[argc-1];

  // Load the solution domain.
  PiecewiseZero dom(inFile, false);
  const int n = dom.length();

  double tau = 2.0/(beta*spring*diffuse);
  double s = 1.0 - exp(-4.0*tim/tau);
  double std = sqrt(s/(beta*spring));
  double denom = sqrt(2.0*M_PI)*std;
  double mu = x0*exp(-2.0*tim/tau);

  printf("diffuse %g\n", diffuse);
  printf("spring %g\n", spring);
  printf("time %g\n", tim);
  printf("beta %g\n", beta);
  printf("x0 %g\n\n", x0);
  printf("tau %g\n", tau);
  printf("s %g\n", s);
  printf("mu %g\n", mu);
  printf("std %g\n", std);

  for (int i = 0; i < n; i++) {
    double x = dom.getNodePos(i);
    double dx = x - mu;
    double p = exp(-0.5*dx*dx/(std*std))/denom;
    dom.set(i, p);
  }
  
  dom.write(outFile);
  return 0;
}

