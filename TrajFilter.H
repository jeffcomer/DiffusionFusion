///////////////////////////////////////////////////////////////////////
// Author: Jeff Comer <jeffcomer at gmail>
// An object for filtering trajectory data.

#ifndef TRAJFILTER_H
#define TRAJFILTER_H

class TrajFilter {
private:
  int cap;
  int minNum;
  int* minVar;
  double* minVal;
  int maxNum;
  int* maxVar;
  double* maxVal;
  int scaleNum;
  int* scaleVar;
  double* scaleVal;
  int periodicNum;
  int* periodicVar;
  double* periodicMin;
  double* periodicMax;
public:
  TrajFilter(int cap0) : cap(cap0), minNum(0), maxNum(0), scaleNum(0), periodicNum(0) {
    minVar = new int[cap];
    minVal = new double[cap];
    maxVar = new int[cap];
    maxVal = new double[cap];
    scaleVar = new int[cap];
    scaleVal = new double[cap];
    periodicVar = new int[cap];
    periodicMin = new double[cap];
    periodicMax = new double[cap];
  }
  ~TrajFilter() {
    delete[] minVar;
    delete[] minVal;
    delete[] maxVar;
    delete[] maxVal;
    delete[] scaleVar;
    delete[] scaleVal;
    delete[] periodicVar;
    delete[] periodicMin;
    delete[] periodicMax;
  }

  double getMinNum() const { return minNum; }
  double getMaxNum() const { return maxNum; }
  double getScaleNum() const { return scaleNum; }
  double getPeriodicNum() const { return periodicNum; }
  // No bounds checking here.
  double getMinVal(int i) const { return minVal[i]; }
  double getMaxVal(int i) const { return maxVal[i]; }
  double getScaleVal(int i) const { return scaleVal[i]; }
  double getPeriodicMin(int i) const { return periodicMin[i]; }
  double getPeriodicMax(int i) const { return periodicMax[i]; }
  int getMinVar(int i) const { return minVar[i]; }
  int getMaxVar(int i) const { return maxVar[i]; }
  int getScaleVar(int i) const { return scaleVar[i]; }
  int getPeriodicVar(int i) const { return periodicVar[i]; }

  bool addMin(int var, double val) {
    if (minNum >= cap) return false;
    minVar[minNum] = var;
    minVal[minNum] = val;
    minNum++;
    return true;
  }
  bool addMax(int var, double val) {
    if (maxNum >= cap) return false;
    maxVar[maxNum] = var;
    maxVal[maxNum] = val;
    maxNum++;
    return true;
  }
  bool addScale(int var, double val) {
    if (scaleNum >= cap) return false;
    scaleVar[scaleNum] = var;
    scaleVal[scaleNum] = val;
    scaleNum++;
    return true;
  }
  bool addPeriodic(int var, double pMin, double pMax) {
    if (periodicNum >= cap) return false;
    periodicVar[periodicNum] = var;
    periodicMin[periodicNum] = pMin;
    periodicMax[periodicNum] = pMax;
    periodicNum++;
    return true;
  }

  String toString() const {
    char s[STRLEN];
    snprintf(s, STRLEN, "filter.min %d filter.max %d filter.scale %d filter.periodic %d", minNum, maxNum, scaleNum, periodicNum);
    return String(s);
  }
private:
  // Don't permit.
  TrajFilter();
  TrajFilter(const TrajFilter&);
  void operator=(const TrajFilter&);
};

#endif
