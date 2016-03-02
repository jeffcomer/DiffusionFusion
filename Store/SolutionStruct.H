// Author: Jeff Comer <jeffcomer at gmail>
#ifndef SOLUTIONSTRUCT_H
#define SOLUTIONSTRUCT_H

struct SolutionStruct {
public:
  int n;
  double dx, dxSq, x0;
  bool periodic;
  double timestep;
};

#endif
