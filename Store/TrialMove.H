// A struct with all public members for holding a descriptor of a node.
// Author: Jeff Comer <jeffcomer at gmail>
#ifndef TRIALMOVE_H
#define TRIALMOVE_H

struct TrialMove {
  int fieldId;
  int node;
  double trialVal, lastVal;
  TrialMove() {}
  TrialMove(int id, int n, double tv, double lv) :
  fieldId(id), node(n), trialVal(tv), lastVal(lv) {
  }
};
#endif
