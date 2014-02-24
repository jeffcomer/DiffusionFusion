///////////////////////////////////////////////////////////////////////
// Author: Jeff Comer <jeffcomer at gmail>
#ifndef EVENT_H
#define EVENT_H

struct Event {
public: 
  static const int varMax = 8;
  double var[varMax]; // variable value
  double del[varMax]; // difference in variable
  double del0[varMax]; // last difference in variable
  int block; // the block number
  int serial; // where does the event lie in the input file?
  int bias; // link to a bias field

  // Organizational parameters put in an unset state.
  Event() : block(-1), serial(-1), bias(-1) {}
};

#endif
