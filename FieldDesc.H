// A struct with all public members for holding a descriptor of a Field.
// Author: Jeff Comer <jeffcomer at gmail>
#ifndef FIELDDESC_H
#define FIELDDESC_H

struct FieldDesc {
  int id;
  String name, type, trajVar;
  bool hasMin, hasMax;
  double minVal, maxVal;
  double step;
  double initVal;
  int n;
  String file;
  bool periodic;
  IndexList varList;
  IndexList priorList;
  bool outLog;
  bool outInt;
  String errFile;
  String fixedFile;
  IndexList freeNodeList;
  bool global;

  String toString() const {
    String s;
    char cs[STRLEN];
    // Usual stuff.
    snprintf(cs, STRLEN, "field %s type %s step %g", name.cs(), type.cs(), step);
    s.add(cs);
    // Optional stuff.
    if (global) {
      snprintf(cs, STRLEN, " global");
      s.add(cs);
    }
    if (file.length() != 0) {
      snprintf(cs, STRLEN, " file %s", file.cs());
      s.add(cs);
    } else {
      snprintf(cs, STRLEN, " n %d initVal %g", n, initVal);
      s.add(cs);
    }
    if (errFile.length() != 0) {
      snprintf(cs, STRLEN, " errFile %s", errFile.cs());
      s.add(cs);
    }
    if (hasMin) {
      snprintf(cs, STRLEN, " minVal %g", minVal);
      s.add(cs);
    }
    if (hasMax) {
      snprintf(cs, STRLEN, " maxVal %g", maxVal);
      s.add(cs);
    }
    if (fixedFile.length() != 0) {
      snprintf(cs, STRLEN, " fixedFile %s", fixedFile.cs());
      s.add(cs);
    }
    if (periodic) {
      snprintf(cs, STRLEN, " periodic");
      s.add(cs);
    }
    for (int v = 0; v < varList.length(); v++) {
      snprintf(cs, STRLEN, " var %d", varList.get(v));
      s.add(cs);
    }
    return s;
  }
};
#endif
