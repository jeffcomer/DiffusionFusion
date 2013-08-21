///////////////////////////////////////////////////////////////////////
// Interface for a Bayesian prior.
// Author: Jeff Comer <jeffcomer at gmail>
#ifndef PRIOR_H
#define PRIOR_H

#include <cstdio>

class Prior {
  const Field* field; // pointer to field master fieldList
  int type;
  int fieldId;
  Field* refField; // a locally stored reference field int gradDim;
  int gradDim;
  double gradVar;

public:
  Prior(String typeName, const Field* field0, int fieldId0) :
    field(field0), type(-1), fieldId(fieldId0), refField(NULL),  gradDim(0), gradVar(-1.0) {

    type = typeIndex(typeName);
    if (type < 0) {
      fprintf(stderr,"ERROR prior: Unknown prior type %s\n", typeName.cs());
      exit(-1);
    }
  }

  ~Prior() {
    if (refField != NULL) delete refField;
  }

  void setRefField(const Field* ref) {
    if (refField != NULL) delete refField;
    if (!field->spannedBy(ref)) {
      fprintf(stderr,"ERROR prior::setRefField refField does not span field %d\n", fieldId);
      fprintf(stderr,"You are probably using a reference field that does not cover the space of the desired field.\n");
      exit(-1);
    }
    // Map the values to the nodes of our field.
    refField = field->map(ref);
  }

  void setRefError(double err) {
    if (refField != NULL) {
      const int n = refField->length();
      for (int i = 0; i < n; i++) refField->setErr(i, err);
    }
  }

  void writeRefField(String fileName) const {
    if (refField != NULL) {
      refField->write(fileName.cs());
      fileName.add(".err");
      refField->writeErr(fileName.cs());
    }
  }

  int getFieldId() { return fieldId; }
  int getType() { return type; }
  String getTypeName() { return name(type); }
  bool hasRefField() { return refField != NULL; }
  
  void setGradDim(int gradDim0) { gradDim=gradDim0; }
  void setGradStd(double gradStd0) {
    if (gradStd0 >= 0) gradVar=gradStd0*gradStd0;
    else gradVar = -1.0;
  }
  int getGradDim() const { return gradDim; }
  double getGradStd() const {
    if (gradVar < 0.0) return -1.0;
    return sqrt(gradVar);
  }

  static const int typeNum = 3;
  static const int scaleType = 0;
  static const int knownType = 1;
  static const int smoothType = 2;
  static const String scaleName;
  static const String knownName;
  static const String smoothName;

  static String name(int type0) {
    switch(type0) {
    case scaleType:
      return scaleName;
    case knownType:
      return knownName;
    case smoothType:
      return smoothName;
    default:
      return String("UNKNOWN_PRIOR");
    }
  }
  static int typeIndex(const String& name) {
    if (name == scaleName) return scaleType;
    if (name == knownName) return knownType;
    if (name == smoothName) return smoothType;
    return -1;
  }

  double calcCost() {
    double cost = 0.0;
    const int n = field->length();

    switch(type) {

      // The scale invariance prior.
    case scaleType:
#pragma omp parallel for reduction(+:cost)
      for (int i = 0; i < n; i++) cost += log(field->get(i));
      break;

      // The known values prior.
    case knownType:
#pragma omp parallel for reduction(+:cost)
      for (int i = 0; i < n; i++) {
	double err = refField->getErr(i);
	double dv = field->get(i) - refField->get(i);
	
	// We'll use negative errors to represent no restraints.
	if (err >= 0.0)
	  cost += 0.5*dv*dv/(err*err);
      }
      break;

      // The smoothness prior.
    case smoothType:
#pragma omp parallel for reduction(+:cost)
      for (int i = 0; i < n; i++) {
	int i0 = field->prevIndex(i, gradDim);
	//printf("gradDim %d i0 %d\n", gradDim, i0);
	// A negative value means nothing to calculate.
	if (i0 >= 0) {
	  double dx = field->spacing(i, gradDim);
	  double dv = field->get(i) - field->get(i0);
	  cost += 0.5*(dv*dv)/(gradVar*dx*dx);
	}
      }
      break;

    }

    return cost;
  }

  double deltaCost(const TrialMove& trialMove) {
    // There is no change if the prior doesn't apply to the modified field.
    if ( trialMove.fieldId != fieldId ) return 0.0;

    int i, i1, i2, i0;
    double err, refV, lastDv, currDv, delta;

    switch(type) {

      // The scale invariance prior.
    case scaleType:
      return log(trialMove.trialVal) - log(trialMove.lastVal);

      // The known values prior.
    case knownType:
      // Check that the change was to this field.
      i = trialMove.node;

      // The known values prior.
      err = refField->getErr(i);
      refV = refField->get(i);
      lastDv = trialMove.lastVal - refV;
      currDv = trialMove.trialVal - refV;
	
      // We'll use negative errors to represent no restraints.
      if (err >= 0.0)
	return 0.5*(currDv*currDv - lastDv*lastDv)/(err*err);
      else return 0.0;

      // The smoothness prior.
    case smoothType:
      delta = 0.0;

      i1 = trialMove.node;
      i0 = field->prevIndex(i1, gradDim);
      i2 = field->nextIndex(i1, gradDim);

      // A negative value means nothing to calculate.
      if (i0 >= 0) {
	double dx = field->spacing(i0, gradDim);
	double v0 = field->get(i0);
	double lastDv = trialMove.lastVal - v0;
	double currDv = trialMove.trialVal - v0;
	delta += 0.5*( currDv*currDv - lastDv*lastDv )/(gradVar*dx*dx);
      }

      // A negative value means nothing to calculate.
      if (i2 >= 0) {
	double dx = field->spacing(i1, gradDim);
	double v2 = field->get(i2);
	double lastDv = v2 - trialMove.lastVal;
	double currDv = v2 - trialMove.trialVal;
	delta += 0.5*( currDv*currDv - lastDv*lastDv )/(gradVar*dx*dx);
      }

      return delta;
    }

    return 0.0;
  }


};

/// Enumeration of the prior types.
const String Prior::scaleName("scale");
const String Prior::knownName("known");
const String Prior::smoothName("smooth");


#endif
