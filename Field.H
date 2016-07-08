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
// Interface for a field of real numbers.
// Author: Jeff Comer <jeffcomer at gmail>
#ifndef FIELD_H
#define FIELD_H

#include "useful.H"

class Field {
public:
  static const int maxDim = 4;
protected:
  double* v0;
  double* err;
  int n;
  bool periodic;

public:
  Field(const double* in, int n0, bool periodic) : n(n0), periodic(false) {
    v0 = new double[n];
    err = new double[n];
    for (int i = 0; i < n; i++) {
      err[i] = 0.0;
      v0[i] = in[i];
    }
  }

  Field(int n0, bool periodic0) : n(n0), periodic(periodic0) {
    v0 = new double[n];
    err = new double[n];
    for (int i = 0; i < n; i++) {
      err[i] = 0.0;
      v0[i] = 0.0;
    }
  }

 Field(const Field& tab) {
    n = tab.n;
    v0 = new double[n];
    err = new double[n];
    for (int i = 0; i < n; i++) {
      v0[i] = tab.v0[i];
      err[i] = tab.err[i];
    }
    periodic = tab.periodic;
  }

  virtual Field& operator=(const Field& tab) {
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
    periodic = tab.periodic;
    
    return (*this);
  }

  virtual ~Field() {
    delete[] v0;
    delete[] err;
  }

  // How many elements?
  int length() const { return n; }
  // The periodicity.
  bool getPeriodic() const { return periodic; }
  void setPeriodic(bool b) { periodic = b; }
  // Get and set elements.
  bool validIndex(int j) const { return (j >= 0) && (j < n); }
  inline double get(int j) const {
    //if (j < 0 || j >= n) return 0.0;
    return v0[j];
  }
  double operator[](int j) const {
    if (j < 0 || j >= n) return 0.0;
    return v0[j];
  }

  // Just dump the numeric data.
  void dump(FILE* out) const {
    for (int i = 0; i < n; i++) fprintf(out, "%.14g\n", v0[i]);
  }

  // Write the error data.
  bool writeErr(const String& fileName) const {
    // Open the file.
    FILE* out = fopen(fileName.cs(),"w");
    if (out == NULL) {
      printf("Couldn't open file `%s' for writing.\n",fileName.cs());
      return false;
    }

    for (int i = 0; i < n; i++) fprintf(out, "%.14g\n", err[i]);
    fclose(out);
    return true;
  }

  void save(double* saveData) const {
    for (int i = 0; i < n; i++) saveData[i] = v0[i];
  }

  // Allow changes to the error.
  double getErr(int j) const {
    if (j < 0 || j >= n) return 0.0;
    return err[j];
  }
  bool setErr(int j, double ev) {
    if (j < 0 || j >= n) return false;
    err[j] = ev;
    return true;
  }
  
  // Basic wrapping if periodic.
  inline int wrapInt(int i, int ni) const {
    if (periodic) {
      int r = i % ni;
      return r<0 ? r+ni : r;
    } else {
      // Aperiodic.
      if (i < 0) return 0;
      if (i >= ni) return ni-1;
      return i;
    }
  }

  // Static wrapping of floating point values, whether the field is periodic or not.
  inline static double wrapReal(double x, double l) {
    int image = int(floor(x/l));
    return x - image*l;
  }
  inline static double wrapReal(double x, double x0, double x1) {
    double l = x1-x0;
    int image = int(floor((x-x0)/l));
    return x - image*l;
  }
  inline static double wrapRealDiff(double x, double l) {
    int image = int(floor(x/l));
    x -= image*l;
    if (x >= 0.5*l) x -= l;
    return x;
  }
  inline static double wrapRealDiff(double x, double x0, double x1) {
    double l = x1-x0;
    int image = int(floor(x/l));
    x -= image*l;
    if (x >= 0.5*l) x -= l;
    return x;
  }

  // Functions that the derived classes must implement.
  // Neighbors must include node j.
  virtual IndexList neighbors(int j) const = 0;
  // Include the minimal neighbors.
  virtual IndexList neighborsMinimal(int j) const = 0;
  // Set the value of a node. Is virtual so that interpolants can be updated.
  virtual bool set(int j, double v) = 0;
  // Read the field from its file.
  virtual void read(const String& fileName, bool periodic0) = 0; 
  // Write the field to a file.
  virtual bool write(const String& fileName) const = 0; 
  // The nearest node to the position (which may be a vector of any size).
  virtual int nearestNode(const double* pos) const = 0;
  // The nearest node with all coordinates less than those of pos.
  virtual int getNode(const double* pos) const = 0;
  // Get the position of a node along dimension dim.
  virtual double nodePos(int j, int dim) const = 0;
  // The integral of the field.
  virtual bool writeIntegral(const String& fileName, double scale) const = 0;
  // The log of the field.
  virtual bool writeLog(const String& fileName, double scale) const = 0;

  // For the smoothness prior.
  virtual int prevIndex(int j, int dir) const = 0;
  virtual int nextIndex(int j, int dir) const = 0;
  virtual double spacing(int j, int dir) const = 0;
  // For the known prior.
  // Map the values of 'srcVal' to the nodes of this.
  // You'll probably have to delete the field returned by this function.
  virtual Field* map(const Field* srcVal) const = 0;
  // Check that 'srcVal' is at least as big as '*this'.
  virtual bool spannedBy(const Field* srcVal) const = 0;

  // How many dimensions does the field have?
  virtual int dimensions() const = 0;

protected:
  // Let PiecewiseCubic take care of all of the initialization.
  Field() {}
};

#endif
