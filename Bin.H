#ifndef BIN_H
#define BIN_H

class Bin {
public:
  Bin() {
  }

  Bin(int n0, double qa, double qb) : n(n0), dx((qb-qa)/n0), x0(qa) {
  }

  static int wrapInt(int i, int m) {
    int r = i % m;
    return r<0 ? r+m : r;
  }

  double getCellPosition(int i) const {
    return x0 + i*dx;
  }

  double getCellCenter(int i) const {
    return x0 + (i+0.5)*dx;
  }

  int getCellWrap(double x) const {
    return wrapInt(int(floor((x - x0)/dx)), n);
  }

  int getCell(double x) const {
    return int(floor((x - x0)/dx));
  }
  
  int getCellNearest(double x) const {
    return int(floor((x - x0)/dx + 0.5));
  }

  int getCellNearestWrap(double x) const {
    return wrapInt(int(floor((x - x0)/dx + 0.5)), n);
  }

public:
  int n;
  double dx;
  double x0;
};

#endif
