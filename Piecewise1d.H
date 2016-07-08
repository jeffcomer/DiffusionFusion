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
#ifndef PIECEWISE1D_H
#define PIECEWISE1D_H

#include "useful.H"
#include "Field.H"

class Piecewise1d : public Field {
protected:
  double dr, dl;
  double r0, r1;
  double e0;
public:

  void read(const String& fileName, bool periodic0) {
    read(fileName.cs(), periodic0);
  }
  void read(const char* fileName, bool periodic0) {
    FILE* inp = fopen(fileName, "r");
    if (inp == NULL) {
      fprintf(stderr,"ERROR Piecewise1d::read Could not open file `%s'\n", fileName);
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
	fprintf(stderr,"ERROR Piecewise1d::read Invalid tabulated potential file line: %s\n", line);
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


  void import(const String& fileSegment, bool periodic0) {
    // Get the data from a string buffer instead of a file.
    int numLines = fileSegment.tokenCount('\n');
    String* lineList = new String[numLines];
    fileSegment.tokenize(lineList, '\n');
    
    double* r = new double[numLines];
    double* v = new double[numLines];
    double* e = new double[numLines];

    int count = 0;
    for (int l = 0; l < numLines; l++) {
      String line = lineList[l].trim();
      int len = lineList[l].length();
      // Ignore comments.
      if (line[0] == '#') continue;
      if (len < 2) continue;
      
      int numTokens = line.tokenCount();
      if (numTokens != 2 && numTokens != 3) {
	fprintf(stderr,"ERROR Piecewise1d::import Invalid import string line: %s\n", line.cs());
	exit(-1);
      }

      String* tokenList = new String[numTokens];
      line.tokenize(tokenList);
      r[count] = strtod(tokenList[0], NULL);
      v[count] = strtod(tokenList[1], NULL);
      if (numTokens == 3) e[count] = strtod(tokenList[2], NULL);
      else e[count] = 0.0;
      count++;

      delete[] tokenList;
    }

    periodic = periodic0;
    init(r, v, count);
    interpolation();

    // Set the error.
    for (int i = 0; i < count; i++) err[i] = e[i];

    delete[] r;
    delete[] v;
    delete[] e;
  }


  //virtual Piecewise1d& operator=(const Piecewise1d& tab) = 0;

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
  virtual void scale(double s) = 0;

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
    
    // Apply a 1d switch.
    double v = v0[indSwitch];
    double m = (value - v)/(indOff - indSwitch);
    for (int i = indSwitch; i < indOff; i++) v0[i] = m*(i - indSwitch) + v;

    interpolation();
    return true;
  }

  virtual double getGrad(int j) const = 0;
  virtual double computeVal(double d) const = 0;
  virtual double computeGrad(double d) const = 0;

  virtual bool computeValGrad(double d, double& val, double& grad) const = 0;

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

  virtual void zero() = 0;

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
  double getDl() const { return dl; }
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
  virtual double nodePos(int j, int dim) const {
    return r0 + dr*j;
  }
  double getNodeCenter(int j) const {
    return r0 + dr*(j+0.5);
  }
  inline int wrapNode(int i) const {
    int r = i % n;
    return r<0 ? r+n : r;
  }

  double wrap(double x) const {
    return wrapReal(x-r0, dl) + r0;
  }
  double wrapDiff(double x) const {
    return wrapRealDiff(x, dl);
  }

  inline int nearestNode(const double* pos) const {
    return nearestNode(pos[0]);
  }
  int nearestNode(double d) const {
    int near;

    if (periodic) {
      near = wrapNode(int(floor((d - r0)/dr + 0.5)));
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
  
  // Check that 'srcVal' is at least as big as '*this'.
  bool spannedBy(const Field* srcVal) const {
    const Piecewise1d* src = dynamic_cast<const Piecewise1d*>(srcVal);
    if (src == NULL) {
      fprintf(stderr,"ERROR Piecewise1d::spannedBy() requires fields of the same type.\n");
      exit(-1);
    }
    printf("%.15g %.15g %.15g %.15g\n", r0+dr, src->r0, r1-dr, src->r1);
    if (r0 + dr < src->r0) return false;
    if (r1 - dr > src->r1) return false;
    return true;
  }

  virtual int dimensions() const { return 1; }

protected:
  // Integrate over cell j.
  // Validation done in integrate(x0,x1).
  // The interpolant gives an exact solution.
  virtual double integrate(int j) const = 0;
  // Integral from the home position rj to x0.
  // Validation done in integrate(x0,x1).
  virtual double integrate(int j, double x) const = 0;

  virtual void init(const double* dist, const double* pot, int n0) = 0;
  

  // Calculate or recalculate the polynomial coefficients
  // for the cubic spline interpolation.
  // We can do this for a range of nodes rather than all.
  virtual void interpolation() = 0;
  
  virtual void interpolation(int j) = 0;

  IndexList neighbors(int j) const {
    IndexList ret;

    if (periodic) {
      for (int k = -2; k <= 1; k++ )
	ret.add(wrapNode(j+k));
    } else {
      for (int k = -2; k <= 1; k++ ) {
	int i = j+k;
	if (i >= 0 && i < n) ret.add(i);
      }
    }

    return ret;
  }

  IndexList neighborsMinimal(int j) const {
    IndexList ret;

    if (periodic) {
      for (int k = -2; k <= 1; k++ )
	ret.add(wrapNode(j+k));
    } else {
      for (int k = -2; k <= 1; k++ ) {
	int i = j+k;
	if (i >= 0 && i < n) ret.add(i);
      }
    }

    return ret;
  }

  virtual void makeInterpolant(int j) = 0;
};
#endif
