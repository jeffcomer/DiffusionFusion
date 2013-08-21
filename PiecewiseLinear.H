///////////////////////////////////////////////////////////////////////
// Author: Jeff Comer <jeffcomer at gmail>
#ifndef PIECEWISELINEAR_H
#define PIECEWISELINEAR_H

#include "useful.H"
#include "Field.H"

class PiecewiseLinear : public Field {
private:
  double* v1;
  double dr, dl;
  double r0, r1;
  double e0;
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

  void read(const String& fileName, bool periodic0) {
    read(fileName.cs(), periodic0);
  }
  void read(const char* fileName, bool periodic0) {
    FILE* inp = fopen(fileName, "r");
    if (inp == NULL) {
      fprintf(stderr,"ERROR PiecewiseLinear::read Could not open file `%s'\n", fileName);
      exit(-1);
    }
    char line[STRLEN];

    const int numLines = countValueLines(fileName);
    double* r = new double[numLines];
    double* v = new double[numLines];
    double* e = new double[numLines];

    int count = 0;
    while (fgets(line, STRLEN, inp) != NULL) {
      // Ignore comments.
      int len = strlen(line);
      if (line[0] == '#') continue;
      if (len < 2) continue;
      
      String s(line);
      int numTokens = s.tokenCount();
      if (numTokens != 2 && numTokens != 3) {
	fprintf(stderr,"ERROR PiecewiseLinear::read Invalid tabulated potential file line: %s\n", line);
	exit(-1);
      }

      String* tokenList = new String[numTokens];
      s.tokenize(tokenList);
      r[count] = strtod(tokenList[0], NULL);
      v[count] = strtod(tokenList[1], NULL);
      if (numTokens == 3) e[count] = strtod(tokenList[2], NULL);
      else e[count] = 0.0;
      count++;

      delete[] tokenList;
    }
    fclose(inp);

    periodic = periodic0;
    init(r, v, count);
    interpolation();

    // Set the error.
    for (int i = 0; i < count; i++) err[i] = e[i];

    delete[] r;
    delete[] v;
    delete[] e;
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

  // We can scale without recalculating the whole interpolant.
  // This should be slightly faster than reinterpolation.
  void scale(double s) {
    for (int i = 0; i < n; i++) {
      v0[i] *= s;
      v1[i] *= s;
      e0 = v1[n-1] + v0[n-1];
    }
  }

  void truncate(double cutoff) {
    int home = int(floor((cutoff - r0)/dr));
    if (home > n) return;

    double v = v0[home];
    for (int i = home; i < n; i++) v0[i] = v;
    interpolation();
  }

  bool truncate(double switchDist, double cutoff, double value) {
    int indOff = int(floor((cutoff - r0)/dr));
    int indSwitch = int(floor((switchDist - r0)/dr));

    if (indSwitch > n) return false;

    // Set everything after the cutoff to "value".
    for (int i = indOff; i < n; i++) v0[i] = value;
    
    // Apply a linear switch.
    double v = v0[indSwitch];
    double m = (value - v)/(indOff - indSwitch);
    for (int i = indSwitch; i < indOff; i++) v0[i] = m*(i - indSwitch) + v;

    interpolation();
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

  double computeZeroOrder(double d) const {
    int near;

    if (periodic) {
      near = int(floor((wrap(d) - r0)/dr + 0.5));
      if (near >= n) near = 0;
    } else {
      near = int(floor((d - r0)/dr + 0.5));
      //printf("INSIDE near %d\n", near);
      if (near < 0) return v0[0];
      else if (near >= n) return e0;
    }

    return v0[near];
  }

  void reset(const double* pot, bool periodic0) {
    periodic = periodic0;
    for (int i = 0; i < n; i++) v0[i] = pot[i];
    interpolation();
  }

  void zero() {
    for (int i = 0; i < n; i++) {
      v0[i] = 0.0;
      v1[i] = 0.0;
    }
  }

  // Write the position and v0 data.
  bool write(const String& fileName) const {
    // Open the file.
    FILE* out = fopen(fileName.cs(),"w");
    if (out == NULL) {
      printf("Couldn't open file %s.\n",fileName.cs());
      return false;
    }

    // Write each point.
    for (int i = 0; i < n; i++) { 
      double r = r0 + i*dr;
      fprintf(out, "%.15g %.15g %.15g\n", r, v0[i], err[i]);
    }

    fclose(out);
    return true;
  }

  // Write the position and the integral of the data.
  bool writeIntegral(const String& fileName, double scale) const {
    double* integ = new double[n];
    
    double integMin = scale*integrate(r0, getNodeCenter(0));
    for (int i = 0; i < n; i++) {
      integ[i] = scale*integrate(r0, getNodeCenter(i));
      if (integ[i] < integMin) integMin = integ[i];
    }
    // Shift the minimum to zero.
    for (int i = 0; i < n; i++) integ[i] -= integMin;

    // Open the file.
    FILE* out = fopen(fileName.cs(),"w");
    if (out == NULL) {
      printf("Couldn't open file %s.\n",fileName.cs());
      return false;
    }

    // Write each point.
    for (int i = 0; i < n; i++) { 
      double r = r0 + (i+0.5)*dr;
      fprintf(out, "%.14g %.14g\n", r, integ[i]);
    }
    fclose(out);

    delete[] integ;
    return true;
  }

  // Write the position and the log of the data
  bool writeLog(const String& fileName, double scale) const {
    // Open the file.
    FILE* out = fopen(fileName.cs(),"w");
    if (out == NULL) {
      printf("Couldn't open file %s.\n",fileName.cs());
      return false;
    }

    // Write each point.
    for (int i = 0; i < n; i++) { 
      double r = r0 + i*dr;
      fprintf(out, "%.14g %.14g\n", r, scale*log(v0[i]));
    }
    fclose(out);

    return true;
  }


  double integrate(double x0, double x1) const {
    if (x0 > x1) return -integrate(x1, x0);
    if (x0 == x1) return 0.0;

    // Integrate zero for values outside of the potential.
    if (x1 < r0) return 0.0;
    if (x0 > r1) return 0.0;
    if (x0 < r0) x0 = r0;
    if (x1 > r1) x1 = r1;
    
    // Find the indices and positions.
    int j0 = int(floor((x0 - r0)/dr));
    int j1 = int(floor((x1 - r0)/dr));
    // We allow integration up to the point r1!
    if (j1 == n) j1 = n-1;
    
    // The normal formula breaks down when j0 == j1.
    if (j0 == j1) return integrate(j0, x1) - integrate(j0, x0);

    // The first segment.
    double ret = integrate(j0) - integrate(j0, x0);
    
    // The internal nodes are easy.
    for (int i = j0+1; i <= j1-1; i++) ret += integrate(i);

    // The last segment.
    ret += integrate(j1, x1);
  
    return ret;
  }

  // Get and set properties.
  double getDr() const { return dr; }
  double getR0() const { return r0; }
  double getR1() const { return r1; }

  bool set(int j, double v) {
    if (j < 0 || j >= n) return false;
    v0[j] = v;

    // We need to update the interpolant.
    // The interpolant is calculated from i-1, i, i+1, and i+2.
    // So from j-2 to j+1 is affected by changing j.
    interpolation(j);
    return true;
  }
  double getNodePos(int j) const {
    return r0 + dr*j;
  }
  double getNodeCenter(int j) const {
    return r0 + dr*(j+0.5);
  }
  inline int wrapNode(int i) const {
    int r = i % n;
    return r<0 ? r+n : r;
  }

  // Wrap a position.
  double wrap(double x) const {
    int image = int(floor( (x-r0)/dl ));
    return x - image*dl;
  }

  inline int nearestNode(const double* pos) const {
    return nearestNode(pos[0]);
  }
  int nearestNode(double d) const {
    int near;

    if (periodic) {
      near = int(floor((wrap(d) - r0)/dr + 0.5));
      if (near >= n) near= 0;
    } else {
      near = int(floor((d - r0)/dr + 0.5));
      if (near < 0) near = 0;
      else if (near >= n) near = n-1;
    }
    return near;
  }

  inline int getNode(const double* pos) const {
    return getNode(pos[0]);
  }
  int getNode(double d) const {
    int near;

    if (periodic) {
      near = int(floor((wrap(d) - r0)/dr));
      if (near >= n) near = 0;
    } else {
      near = int(floor((d - r0)/dr));
      if (near < 0) near = 0;
      else if (near >= n) near = n-1;
    }
    return near;
  }

  // Return the index of the previous node along a dimension.
  // Return -1 for none.
  int prevIndex(int j, int dim) const {
    // We don't check that dim==0.
    if (periodic) return wrapNode(j-1);
    else {
      // Ignore out-of-bounds or j=0.
      if (j <= 0 || j >= n) return -1;
      return j-1;
    }
  }
  // Return the index of the next node along a dimension.
  // Return -1 for none.
  int nextIndex(int j, int dim) const {
    // We don't check dim==0.
    if (periodic) return wrapNode(j+1);
    else {
      // Ignore out-of-bounds or j=n-1.
      if (j < 0 || j >= n-1) return -1;
      return j+1;
    }
  }
  // Grid spacing.
  double spacing(int j, int dim) const { return dr; }
  
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

  // Check that 'srcVal' is at least as big as '*this'.
  bool spannedBy(const Field* srcVal) const {
    const PiecewiseLinear* src = dynamic_cast<const PiecewiseLinear*>(srcVal);
    if (src == NULL) {
      fprintf(stderr,"ERROR PiecewiseLinear::spannedBy() requires fields of the same type.\n");
      exit(-1);
    }
    printf("%.15g %.15g %.15g %.15g\n", r0+dr, src->r0, r1-dr, src->r1);
    if (r0 + dr < src->r0) return false;
    if (r1 - dr > src->r1) return false;
    return true;
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

  IndexList neighbors(int j) const {
    IndexList ret;

    if (periodic) {
      for (int k = -2; k <= 2; k++ )
	ret.add(wrapNode(j+k));
    } else {
      for (int k = -2; k <= 2; k++ ) {
	int i = j+k;
	if (i >= 0 && i < n) ret.add(i);
      }
    }

    return ret;
  }

  IndexList neighborsWide(int j) const {
    IndexList ret;

    if (periodic) {
      for (int k = -2; k <= 2; k++ )
	ret.add(wrapNode(j+k));
    } else {
      for (int k = -2; k <= 2; k++ ) {
	int i = j+k;
	if (i >= 0 && i < n) ret.add(i);
      }
    }

    return ret;
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
