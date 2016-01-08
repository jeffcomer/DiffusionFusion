///////////////////////////////////////////////////////////////////////
// Author: Jeff Comer <jeffcomer at gmail>
#ifndef PIECEWISELINEAR_H
#define PIECEWISELINEAR_H

#include "useful.H"
#include "Piecewise1d.H"

class PiecewiseLinear : public Piecewise1d {
private:
  double* v1;
public:
  PiecewiseLinear(const char* fileName, bool periodic0) {
    v0 = NULL;
    read(fileName, periodic0);
  }

  PiecewiseLinear(const PiecewiseLinear& tab) {
    n = tab.n;
    dr = tab.dr;
    dl = tab.dl;
    r0 = tab.r0;
    r1 = tab.r1;

    v0 = new double[n];
    v1 = new double[n];
    err = new double[n];
    for (int i = 0; i < n; i++) {
      v0[i] = tab.v0[i];
      v1[i] = tab.v1[i];
      err[i] = tab.err[i];
    }
    e0 = tab.e0;
    periodic = tab.periodic;
  }
  
  PiecewiseLinear(const Piecewise1d& tab) {
    n = tab.length();
    dr = tab.getDr();
    dl = tab.getDl();
    r0 = tab.getR0();
    r1 = tab.getR1();
    periodic = tab.getPeriodic();
    
    v0 = new double[n];
    v1 = new double[n];
    err = new double[n];
    for (int i = 0; i < n; i++) {
      v0[i] = tab.get(i);
      err[i] = tab.getErr(i);
    }
    interpolation();
  }

  PiecewiseLinear(const double* dist, const double* pot, int n0, bool periodic0) {
    v0 = NULL;
    periodic = periodic0;
    init(dist, pot, n0);
    interpolation();
  }

  PiecewiseLinear(double x0, double x1, double pot0, int n0, bool periodic0) {
    n = abs(n0);
    r0 = x0;
    r1 = x1;
    dl = r1-r0;
    dr = dl/n;
    periodic = periodic0;
    
    v0 = new double[n];
    v1 = new double[n];
    err = new double[n];
    for (int i = 0; i < n; i++) {
      v0[i] = pot0;
      err[i] = 0.0;
    }
    interpolation();
  }

  ~PiecewiseLinear() {
    delete[] v1;
  }

  PiecewiseLinear& operator=(const PiecewiseLinear& tab) {
    // Resize the arrays if necessary.
    if (n != tab.n) {
      delete[] v0;
      delete[] v1;
      delete[] err;
      n = tab.n;
      v0 = new double[n];
      v1 = new double[n];
      err = new double[n];
    }

    // Copy the values.
    for (int i = 0; i < n; i++) {
      v0[i] = tab.v0[i];
      v1[i] = tab.v1[i];
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

  // Get the gradient value.
  double getGrad(int j) const {
    if (j < 0 || j >= n) return 0.0;
    return v1[j]/dr;
  }

  // We can scale without recalculating the whole interpolant.
  // This should be slightly faster than reinterpolation.
  void scale(double s) {
    for (int i = 0; i < n; i++) {
      v0[i] *= s;
      v1[i] *= s;
    }
    e0 = v1[n-1] + v0[n-1];
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
   
    // Interpolate.
    return v1[home]*w + v0[home];
  }

  double computeGrad(double d) const {
    int home;

    if (periodic) {
      home = int(floor((wrap(d) - r0)/dr));
    } else {
      home = int(floor((d - r0)/dr));
      //printf("INSIDE home %d\n", home);
      if (home < 0) return v0[0];
      else if (home >= n) return e0;
    }
  
    //double homeR = home*dr + r0;
    //double w = (d - homeR)/dr;
   
    // Zero-order.
    return v1[home];
  }

  bool computeValGrad(double d, double& val, double& grad) const {
    val =  computeVal(d);
    grad =  computeGrad(d);
    return true;
  }

  void zero() {
    for (int i = 0; i < n; i++) {
      v0[i] = 0.0;
      v1[i] = 0.0;
    }
  }

  // Map. Don't forget to delete the returned value.
  // Map the values of 'srcVal' to the nodes of this.
  Field* map(const Field* srcVal) const {
    const PiecewiseLinear* src = dynamic_cast<const PiecewiseLinear*>(srcVal);
    if (src == NULL) {
      fprintf(stderr,"ERROR PiecewiseLinear::map() requires fields of the same type.\n");
      exit(-1);
    }
    
    PiecewiseLinear* ret = new PiecewiseLinear(*this); 
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
    return dr*(0.5*v1[j] + v0[j]);
  }
  // Integral from the home position rj to x0.
  // Validation done in integrate(x0,x1).
  double integrate(int j, double x) const {
    double rj = r0 + dr*j;
    double w = (x - rj)/dr;
    double w2 = w*w;

    return dr*(0.5*v1[j]*w2 + v0[j]*w);
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
      delete[] err;
    }
    
    v0 = new double[n];
    v1 = new double[n];
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
    e0 = v1[n-1] + v0[n-1];
  }
  
  void interpolation(int j) { 
    // The interpolant is calculated from i, i+1.
    // So from j-1 to j are affected by changing j.

    if (periodic) {
      for (int k = -2; k <= 2; k++ ) makeInterpolant(wrapNode(j+k));
    } else {
      for (int k = -2; k <= 2; k++ ) {
	int i = j+k;
	if (i >= 0 && i < n) makeInterpolant(i);
      }
    }
    e0 = v1[n-1] + v0[n-1];
  }

  void makeInterpolant(int j) {
    // We assume that j is a valid node.
    int i1 = j;
    int i2 = j+1;
    
    // Use periodic or constant boundary conditions.
    if (periodic) {
      if (i2 >= n) i2 = 0;
    } else {
      if (i2 >= n) i2 = n-1;
    }

    // Calculate the interpolant parameters.
    v1[j] = v0[i2] - v0[i1];
  }
};
#endif
