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
// Interface for a Bayesian prior.
// Author: Jeff Comer <jeffcomer at gmail>
#ifndef PRIOR_H
#define PRIOR_H

#include <cstdio>
#include "TrialMove.H"

class Prior {
  const Field* field; // pointer to field master fieldList
  int type;
  int fieldId;
  Field* refField; // a locally stored reference field int gradDim;
  int gradDim;
  double gradVar;
  int startIndex, endIndex;
  const Field* coupleField;
  int coupleFieldId;
  double coupleVar;

public:
  Prior(String typeName, const Field* field0, int fieldId0) :
    field(field0), type(-1), fieldId(fieldId0), refField(NULL),  gradDim(0), gradVar(-1.0),
    coupleField(NULL), coupleVar(-1.0) {

    startIndex = 0;
    endIndex = field->length() - 1;

    type = typeIndex(typeName);
    if (type < 0) {
      fprintf(stderr,"ERROR prior: Unknown prior type %s\n", typeName.cs());
      exit(-1);
    }
  }

  ~Prior() {
    if (refField != NULL) delete refField;
    if (coupleField != NULL) delete coupleField;
  }


  void setCoupleField(const Field* couple, int coupleFieldId0) {
    if (coupleField != NULL) delete coupleField;

    coupleField = couple;
    if (coupleField->length() != field->length()) {
      fprintf(stderr,"ERROR prior::setCoupleField coupleField does not have the same number of nodes as field %d\n", fieldId);
      exit(-1);
    }
    coupleFieldId = coupleFieldId0;
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

  bool setStartIndex(int ind) {
    if (ind >= 0 && ind < field->length()) {
      startIndex = ind;
      return true;
    }
    return false;
  }
  bool setEndIndex(int ind) {
    if (ind >= 0 && ind < field->length()) {
      endIndex = ind;
      return true;
    }
    return false;
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
  bool hasCoupleField() { return coupleField != NULL; }
  
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

  void setCoupleStd(double coupleStd) {
    if (coupleStd >= 0) coupleVar = coupleStd*coupleStd;
  }
  double getCoupleStd() const {
    if (coupleVar < 0.0) return -1.0;
    return sqrt(coupleVar);
  }

  static const int typeNum = 4;
  static const int scaleType = 0;
  static const int knownType = 1;
  static const int smoothType = 2;
  static const int coupleType = 3;
  static const String scaleName;
  static const String knownName;
  static const String smoothName;
  static const String coupleName;

  static String name(int type0) {
    switch(type0) {
    case scaleType:
      return scaleName;
    case knownType:
      return knownName;
    case smoothType:
      return smoothName;
    case coupleType:
      return coupleName;
    default:
      return String("UNKNOWN_PRIOR");
    }
  }
  static int typeIndex(const String& name) {
    if (name == scaleName) return scaleType;
    if (name == knownName) return knownType;
    if (name == smoothName) return smoothType;
    if (name == coupleName) return coupleType;
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
      if (refField == NULL) {
	fprintf(stderr, "ERROR Prior known: reference field has not been set.\n");
	  exit(-1);
      }
#pragma omp parallel for reduction(+:cost)
      for (int i = startIndex; i <= endIndex; i++) {
	double err = refField->getErr(i);
	double dv = field->get(i) - refField->get(i);
	
	// We'll use negative errors to represent no restraints.
	if (err >= 0.0)
	  cost += 0.5*dv*dv/(err*err);
      }
      break;

      // The smoothness prior.
    case smoothType:
      // A pair is considered valid for the smoothness prior when the node with
      // the larger index is consistent with startIndex <= index <= endIndex.
#pragma omp parallel for reduction(+:cost)
      for (int i = startIndex; i <= endIndex; i++) {
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

     // The smoothness prior.
    case coupleType:
      if (coupleField == NULL) {
	fprintf(stderr, "ERROR Prior couple: field to couple to has not been set.\n");
	exit(-1);
      }

      // A pair is considered valid for the smoothness prior when the node with
      // the larger index is consistent with startIndex <= index <= endIndex.
#pragma omp parallel for reduction(+:cost)
      for (int i = startIndex; i <= endIndex; i++) {
	double dv = field->get(i) - coupleField->get(i);
	cost += 0.5*(dv*dv)/coupleVar;
      }
      break;

    }

    return cost;
  }

  double deltaCost(const TrialMove& trialMove) {
    int i, i1, i2, i0;
    double err, refV, lastDv, currDv, delta;

    switch(type) {


    case scaleType:
      // The scale invariance prior.
      // There is no change if the prior doesn't apply to the modified field.
      if ( trialMove.fieldId != fieldId ) return 0.0;
      return log(trialMove.trialVal) - log(trialMove.lastVal);


    case knownType:
      // The known values prior.
      // There is no change if the prior doesn't apply to the modified field.
      if ( trialMove.fieldId != fieldId ) return 0.0;

      i = trialMove.node;

      // The known values prior.
      err = refField->getErr(i);
      refV = refField->get(i);
      lastDv = trialMove.lastVal - refV;
      currDv = trialMove.trialVal - refV;
	
      // We'll use negative errors to represent no restraints.
      if (err >= 0.0 && i >= startIndex && i <= endIndex )
	return 0.5*(currDv*currDv - lastDv*lastDv)/(err*err);
      else return 0.0;

    case smoothType:
      // The smoothness prior.
      // There is no change if the prior doesn't apply to the modified field.
      if ( trialMove.fieldId != fieldId ) return 0.0;
      delta = 0.0;

      i1 = trialMove.node;
      i0 = field->prevIndex(i1, gradDim);
      i2 = field->nextIndex(i1, gradDim);
      // A pair is considered valid for the smoothness prior when the node with
      // the larger index is consistent with startIndex <= index <= endIndex.
      
      // A negative value means nothing to calculate.
      if (i0 >= 0 && i1 >= startIndex && i1 <= endIndex ) {
	double dx = field->spacing(i0, gradDim);
	double v0 = field->get(i0);
	double lastDv = trialMove.lastVal - v0;
	double currDv = trialMove.trialVal - v0;
	delta += 0.5*( currDv*currDv - lastDv*lastDv )/(gradVar*dx*dx);
      }

      // A negative value means nothing to calculate.
      if (i2 >= 0 && i2 >= startIndex && i2 <= endIndex ) {
	double dx = field->spacing(i1, gradDim);
	double v2 = field->get(i2);
	double lastDv = v2 - trialMove.lastVal;
	double currDv = v2 - trialMove.trialVal;
	delta += 0.5*( currDv*currDv - lastDv*lastDv )/(gradVar*dx*dx);
      }

      return delta;

    case coupleType:
      // The coupling prior.
      // There is no change if the prior doesn't apply to the modified field.
      // We must check both the prior's field and field it's COUPLED to!
      if ( trialMove.fieldId != fieldId && trialMove.fieldId != coupleFieldId ) return 0.0;

      i = trialMove.node;
      if (i < startIndex || i > endIndex ) return 0.0;

      // Coupling prior.
      double v;
      if (trialMove.fieldId == fieldId) {
	// Move was in the primary field.
	v = coupleField->get(i);
      } else {
	// Move was in coupleField.
	v = field->get(i);
      }
      lastDv = trialMove.lastVal - v;
      currDv = trialMove.trialVal - v;

      return 0.5*(currDv*currDv - lastDv*lastDv)/coupleVar;
    }

    return 0.0;
  }


};

/// Enumeration of the prior types.
const String Prior::scaleName("scale");
const String Prior::knownName("known");
const String Prior::smoothName("smooth");
const String Prior::coupleName("couple");

#endif
