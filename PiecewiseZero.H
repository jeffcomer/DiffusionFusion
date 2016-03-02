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
#ifndef PIECEWISEZERO_H
#define PIECEWISEZERO_H

#include "useful.H"
#include "Field.H"
#include "Piecewise1d.H"

class PiecewiseZero : public Piecewise1d {
public:
  PiecewiseZero(const char* fileName, bool periodic0) {
    v0 = NULL;
    read(fileName, periodic0);
  }

  PiecewiseZero(const PiecewiseZero& tab) {
    n = tab.n;
    dr = tab.dr;
    dl = tab.dl;
    r0 = tab.r0;
    r1 = tab.r1;

    v0 = new double[n];
    err = new double[n];
    for (int i = 0; i < n; i++) {
      v0[i] = tab.v0[i];
      err[i] = tab.err[i];
    }
    e0 = tab.e0;
    periodic = tab.periodic;
  }
  
  PiecewiseZero(const double* dist, const double* pot, int n0, bool periodic0) {
    v0 = NULL;
    periodic = periodic0;
    init(dist, pot, n0);
    e0 = v0[n-1];
  }

  PiecewiseZero(const double* pot, int n0, double x00, double dx0, bool periodic0) {
    n = abs(n0);
    dr = fabs(dx0);
    r0 = x00;
    r1 = r0 + n*dr;
    dl = r1 - r0;

    v0 = new double[n];
    err = new double[n];
    for (int i = 0; i < n; i++) {
      v0[i] = pot[i];
      err[i] = 0.0;
    }
    e0 = v0[n-1];
  }

  PiecewiseZero(double x0, double x1, double pot0, int n0, bool periodic0) {
    n = abs(n0);
    r0 = x0;
    r1 = x1;
    dl = r1-r0;
    dr = dl/n;
    periodic = periodic0;
    
    v0 = new double[n];
    err = new double[n];
    for (int i = 0; i < n; i++) {
      v0[i] = pot0;
      err[i] = 0.0;
    }
    e0 = v0[n-1];
  }

  ~PiecewiseZero() {
  }

  PiecewiseZero& operator=(const PiecewiseZero& tab) {
    // Resize the arrays if necessary.
    if (n != tab.n) {
      delete[] v0;
      delete[] err;
      n = tab.n;
      v0 = new double[n];
      err = new double[n];
    }

    // Copy the values.
    for (int i = 0; i < n; i++) {
      v0[i] = tab.v0[i];
      err[i] = tab.err[i];
    }
    dr = tab.dr;
    dl = tab.dl;
    r0 = tab.r0;
    r1 = tab.r1;
    e0 = tab.e0;
    periodic = tab.periodic;
    
    return (*this);
  }

  static int countValueLines(const char* fileName) {
    FILE* inp = fopen(fileName, "r");
    char line[STRLEN];
    int count = 0;

    while (fgets(line, STRLEN, inp) != NULL) {
      // Ignore comments.
      int len = strlen(line);
      if (line[0] == '#') continue;
      if (len < 2) continue;
      
      count++;
    }
    fclose(inp);

    return count;
  }

  double getGrad(int j) const {
    if (j < 0 || j >= n) return 0.0;
    if (j == 0) return (v0[1]-v0[0])/dr;
    if (j == n-1) return (v0[n-1]-v0[n-2])/dr;
    return (v0[j+1]-v0[j-1])/(2.0*dr);
  }

  void scale(double s) {
    for (int i = 0; i < n; i++) {
      v0[i] *= s;
    }
  }

  double computeVal(double d) const {
    int home;

    if (periodic) {
      home = int(floor((wrap(d) - r0)/dr));
    } else {
      home = int(floor((d - r0)/dr));
      //printf("INSIDE home %d\n", home);
      if (home < 0) return v0[0];
      else if (home >= n) return v0[n-1];
    }

    return v0[home];
  }

  double computeGrad(double d) const {
    return 0.0;
  }

  bool computeValGrad(double d, double& val, double& grad) const {
    val =  computeVal(d);
    grad = 0.0;
    return true;
  }

  void zero() {
    for (int i = 0; i < n; i++) {
      v0[i] = 0.0;
    }
  }

  // Map. Don't forget to delete the returned value.
  // Map the values of 'srcVal' to the nodes of this.
  Field* map(const Field* srcVal) const {
    const PiecewiseZero* src = dynamic_cast<const PiecewiseZero*>(srcVal);
    if (src == NULL) {
      fprintf(stderr,"ERROR PiecewiseZero::map() requires fields of the same type.\n");
      exit(-1);
    }
    
    PiecewiseZero* ret = new PiecewiseZero(*this); 
    for (int i = 0; i < n; i++) {
      double x = getNodePos(i);
      double v = src->computeVal(x);
      int near = src->nearestNode(x);
      ret->v0[i] = v;
      ret->err[i] = src->err[near];
    }

    return ret;
  }


private:
  // Integrate over cell j.
  // Validation done in integrate(x0,x1).
  // The interpolant gives an exact solution.
  double integrate(int j) const {
    return dr*v0[j];
  }
  // Integral from the home position rj to x0.
  // Validation done in integrate(x0,x1).
  double integrate(int j, double x) const {
    double rj = r0 + dr*j;
    double w = (x - rj)/dr;

    return dr*v0[j]*w;
  }

  void init(const double* dist, const double* pot, int n0) {
    n = abs(n0);
    dr = abs(dist[1]-dist[0]);
    r0 = dist[0]<=dist[1] ? dist[0] : dist[1];
    r1 = r0 + n*dr;
    dl = r1 - r0;

    if (v0 != NULL) {
      delete[] v0;
      delete[] err;
    }
    
    v0 = new double[n];
    err = new double[n];
    for (int i = 0; i < n; i++) {
      v0[i] = pot[i];
      err[i] = 0.0;
    }
  }
  
  // Interpolation functions do nothing since this is zero-order.
  void interpolation() {}
  void interpolation(int j) {}
  void makeInterpolant(int j) {}
};
#endif
