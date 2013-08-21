///////////////////////////////////////////////////////////////////////
// Author: Jeff Comer <jeffcomer at gmail>
#ifndef EVENT_H
#define EVENT_H

class Event {
public: 
  static const int varMax = 8;
  double var[varMax]; // variable value
  double del[varMax]; // difference in variable
  double del0[varMax]; // last difference in variable
};

#endif
