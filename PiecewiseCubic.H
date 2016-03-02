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
#ifndef PIECEWISECUBIC_H
#define PIECEWISECUBIC_H

#include "useful.H"
#include "Piecewise1d.H"

class PiecewiseCubic : public Piecewise1d {
private:
  double* v1;
  double* v2;
  double* v3;
public:
  PiecewiseCubic(const char* fileName, bool periodic0) {
    v0 = NULL;
    read(fileName, periodic0);
  }

  PiecewiseCubic(bool periodic0, const String& fileSegment) {
    v0 = NULL;
    import(fileSegment,periodic0);
  }

  PiecewiseCubic(const PiecewiseCubic& tab) {
    n = tab.n;
    dr = tab.dr;
    dl = tab.dl;
    r0 = tab.r0;
    r1 = tab.r1;

    v0 = new double[n];
    v1 = new double[n];
    v2 = new double[n];
    v3 = new double[n];
    err = new double[n];
    for (int i = 0; i < n; i++) {
      v0[i] = tab.v0[i];
      v1[i] = tab.v1[i];
      v2[i] = tab.v2[i];
      v3[i] = tab.v3[i];
      err[i] = tab.err[i];
    }
    e0 = tab.e0;
    periodic = tab.periodic;
  }
  
  PiecewiseCubic(const Piecewise1d& tab) {
    n = tab.length();
    dr = tab.getDr();
    dl = tab.getDl();
    r0 = tab.getR0();
    r1 = tab.getR1();
    periodic = tab.getPeriodic();
    
    v0 = new double[n];
    v1 = new double[n];
    v2 = new double[n];
    v3 = new double[n];
    err = new double[n];
    for (int i = 0; i < n; i++) {
      v0[i] = tab.get(i);
      err[i] = tab.getErr(i);
    }
    interpolation();
  }
  
  PiecewiseCubic(const double* dist, const double* pot, int n0, bool periodic0) {
    v0 = NULL;
    periodic = periodic0;
    init(dist, pot, n0);
    interpolation();
  }

  PiecewiseCubic(double x0, double x1, const double* pot, int n0, bool periodic0) {
    n = abs(n0);
    r0 = x0;
    r1 = x1;
    dl = r1-r0;
    dr = dl/n;
    periodic = periodic0;
    
    v0 = new double[n];
    v1 = new double[n];
    v2 = new double[n];
    v3 = new double[n];
    err = new double[n];
    for (int i = 0; i < n; i++) {
      v0[i] = pot[i];
      err[i] = 0.0;
    }
    interpolation();
  }

  PiecewiseCubic(double x0, double x1, double pot0, int n0, bool periodic0) {
    n = abs(n0);
    r0 = x0;
    r1 = x1;
    dl = r1-r0;
    dr = dl/n;
    periodic = periodic0;
    
    v0 = new double[n];
    v1 = new double[n];
    v2 = new double[n];
    v3 = new double[n];
    err = new double[n];
    for (int i = 0; i < n; i++) {
      v0[i] = pot0;
      err[i] = 0.0;
    }
    interpolation();
  }

  ~PiecewiseCubic() {
    delete[] v1;
    delete[] v2;
    delete[] v3;
  }


  PiecewiseCubic& operator=(const PiecewiseCubic& tab) {
    // Resize the arrays if necessary.
    if (n != tab.n) {
      delete[] v0;
      delete[] v1;
      delete[] v2;
      delete[] v3;
      delete[] err;
      n = tab.n;
      v0 = new double[n];
      v1 = new double[n];
      v2 = new double[n];
      v3 = new double[n];
      err = new double[n];
    }

    // Copy the values.
    for (int i = 0; i < n; i++) {
      v0[i] = tab.v0[i];
      v1[i] = tab.v1[i];
      v2[i] = tab.v2[i];
      v3[i] = tab.v3[i];
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


  // We can scale without recalculating the whole interpolant.
  // This should be slightly faster than reinterpolation.
  void scale(double s) {
    for (int i = 0; i < n; i++) {
      v0[i] *= s;
      v1[i] *= s;
      v2[i] *= s;
      v3[i] *= s;
    }
    e0 = v3[n-1] + v2[n-1] + v1[n-1] + v0[n-1];
  }

  // We can shift without recalculating the whole interpolant.
  // This should be slightly faster than reinterpolation.
  void shift(double s) {
    for (int i = 0; i < n; i++) {
      v0[i] += s;
    }
    e0 = v3[n-1] + v2[n-1] + v1[n-1] + v0[n-1];
  }

  // Get the gradient value.
  inline double getGrad(int j) const {
    if (j < 0 || j >= n) return 0.0;
    if (j == 0) return (v0[1]-v0[0])/dr;
    if (j == n-1) return (v0[n-1]-v0[n-2])/dr;
    return v1[j]/dr;
  }

 double computeGrad(double d) const {
    int home;

    if (periodic) {
      d = wrap(d);
      home = int(floor((d - r0)/dr));
    } else {
      home = int(floor((d - r0)/dr));

      if (home < 0) return 0.0;
      if (home >= n) return 0.0;
    }

    double homeR = home*dr + r0;
    double w = (d - homeR)/dr;
    double w2 = w*w;
   
    // Interpolate.
    return (3.0*v3[home]*w2 + 2.0*v2[home]*w + v1[home])/dr;
  }

 double computeGradGrad(double d) const {
    int home;

    if (periodic) {
      d = wrap(d);
      home = int(floor((d - r0)/dr));
    } else {
      home = int(floor((d - r0)/dr));

      if (home < 0) return 0.0;
      if (home >= n) return 0.0;
    }

    double homeR = home*dr + r0;
    double w = (d - homeR)/dr;
   
    // Interpolate.
    return (6.0*v3[home]*w + 2.0*v2[home])/(dr*dr);
  }

  bool computeValGrad(double d, double& val, double& grad) const {
    int home;

    if (periodic) {
      d = wrap(d);
      home = int(floor((d - r0)/dr));
    } else {
      home = int(floor((d - r0)/dr));

      if (home < 0) {
	val = v0[0];
	grad = 0.0;
	return false;
      } else if (home >= n) {
	val = e0;
	grad = 0.0;
	return false;
      }
    }

    double homeR = home*dr + r0;
    double w = (d - homeR)/dr;
    double w2 = w*w;
   
    // Interpolate.
    val = v3[home]*w*w2 + v2[home]*w2 + v1[home]*w + v0[home];
    grad = (3.0*v3[home]*w2 + 2.0*v2[home]*w + v1[home])/dr;
    return true;
  }

  double computeVal(double d) const {
    int home;

    if (periodic) {
      home = int(floor((wrap(d) - r0)/dr));
    } else {
      home = int(floor((d - r0)/dr));
      //printf("INSIDE home %d\n", home);
      if (home < 0) return v0[0];
      else if (home >= n) return e0;
    }
  
    double homeR = home*dr + r0;
    double w = (d - homeR)/dr;
    double w2 = w*w;
   
    // Interpolate.
    return v3[home]*w*w2 + v2[home]*w2 + v1[home]*w + v0[home];
  }

  void zero() {
    for (int i = 0; i < n; i++) {
      v0[i] = 0.0;
      v1[i] = 0.0;
      v2[i] = 0.0;
      v3[i] = 0.0;
    }
  }

  // Map. Don't forget to delete the returned value.
  // Map the values of 'srcVal' to the nodes of this.
  Field* map(const Field* srcVal) const {
    const PiecewiseCubic* src = dynamic_cast<const PiecewiseCubic*>(srcVal);
    if (src == NULL) {
      fprintf(stderr,"ERROR PiecewiseCubic::map() requires fields of the same type.\n");
      exit(-1);
    }
    
    PiecewiseCubic* ret = new PiecewiseCubic(*this); 
    for (int i = 0; i < n; i++) {
      double x = getNodePos(i);
      double v = src->computeVal(x);
      int near = src->nearestNode(x);
      ret->v0[i] = v;
      ret->err[i] = src->err[near];
    }
    ret->interpolation();

    return ret;
  }


private:
  // Integrate over cell j.
  // Validation done in integrate(x0,x1).
  // The interpolant gives an exact solution.
  double integrate(int j) const {
    //if (j < 0 || j >= n) return 0.0;
    return dr*(0.25*v3[j] + v2[j]/3.0 + 0.5*v1[j] + v0[j]);
  }
  // Integral from the home position rj to x0.
  // Validation done in integrate(x0,x1).
  double integrate(int j, double x) const {
    double rj = r0 + dr*j;
    double w = (x - rj)/dr;
    double w2 = w*w;

    return dr*(0.25*v3[j]*w2*w2 + v2[j]/3.0*w2*w + 0.5*v1[j]*w2 + v0[j]*w);
  }

  void init(const double* dist, const double* pot, int n0) {
    n = abs(n0);
    dr = abs(dist[1]-dist[0]);
    r0 = dist[0]<=dist[1] ? dist[0] : dist[1];
    r1 = r0 + n*dr;
    dl = r1 - r0;

    if (v0 != NULL) {
      delete[] v0;
      delete[] v1;
      delete[] v2;
      delete[] v3;
      delete[] err;
    }
    
    v0 = new double[n];
    v1 = new double[n];
    v2 = new double[n];
    v3 = new double[n];
    err = new double[n];
    for (int i = 0; i < n; i++) {
      v0[i] = pot[i];
      err[i] = 0.0;
    }
  }
  

  // Calculate or recalculate the polynomial coefficients
  // for the cubic spline interpolation.
  // We can do this for a range of nodes rather than all.
  void interpolation() {
    for (int j = 0; j < n; j++) makeInterpolant(j);
    e0 = v3[n-1] + v2[n-1] + v1[n-1] + v0[n-1];
  }
  
  void interpolation(int j) { 
    // The interpolant is calculated from i-1, i, i+1, and i+2.
    // So from j-2 to j+1 are affected by changing j.

    if (periodic) {
      for (int k = -2; k <= 1; k++ ) makeInterpolant(wrapNode(j+k));
    } else {
      for (int k = -2; k <= 1; k++ ) {
	int i = j+k;
	if (i >= 0 && i < n) makeInterpolant(i);
      }
    }
    e0 = v3[n-1] + v2[n-1] + v1[n-1] + v0[n-1];
  }

  void makeInterpolant(int j) {
    // We assume that j is a valid node.
    int i0 = j - 1;
    int i1 = j;
    int i2 = j + 1;
    int i3 = j + 2;
    
    // Use periodic or constant boundary conditions.
    if (periodic) {
      if (i0 < 0) i0 += n;
      if (i2 >= n) i2 -= n;
      if (i3 >= n) i3 -= n;
    } else {
      if (i0 < 0) i0 = 0;
      if (i2 >= n) i2 = n-1;
      if (i3 >= n) i3 = n-1;
    }

    // Calculate the interpolant parameters.
    v3[j] = 0.5*(-v0[i0] + 3*v0[i1] - 3*v0[i2] + v0[i3]);
    v2[j] = 0.5*(2*v0[i0] - 5*v0[i1] + 4*v0[i2] - v0[i3]);
    v1[j] = 0.5*(-v0[i0] + v0[i2]);
  }
};
#endif
