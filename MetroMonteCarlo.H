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
// Metropolis Monte Carlo over fields
// Author: Jeff Comer <jeffcomer at gmail>
#ifndef METROMONTECARLO_H
#define METROMONTECARLO_H

#include "Field.H"
#include "FieldDesc.H"
#include "TrialMove.H"
#include "RandomGsl.H"

class MetroMonteCarlo {
private: 
  // For the fields on which we are applying Monte Carlo.
  int fieldNum;
  // An array of pointers to the fields on which we are applying Monte Carlo.
  Field** mcField;
  // An array of pointers to their descriptors.
  const FieldDesc** mcFieldDesc;
  // Random number generator
  Random* rando;
  IndexList mcLocalFields;

  // Trial move variables.
  int trialField, trialNode;
  double trialLastVal, trialVal;

  int stepsPerCycle;

public:
  MetroMonteCarlo(Field** ppf, const FieldDesc* pfd, const IndexList& selection, Random* rando0) {
    fieldNum = selection.length();
    if (fieldNum == 0) {
      fprintf(stderr, "ERROR MetroMonteCarlo must act on at least one field.\n");
      exit(-1);
    }

    // Choose the selection of fields to undergo Monte Carlo.
    mcField = new Field*[fieldNum];
    mcFieldDesc = new const FieldDesc*[fieldNum];
    for (int i = 0; i < fieldNum; i++) {
      mcField[i] = ppf[selection.get(i)];
      mcFieldDesc[i] = &pfd[selection.get(i)];
      if (!mcFieldDesc[i]->global) mcLocalFields.add(i);
    }

    // Get the number of steps per cycle (equal to the number of field nodes).
    stepsPerCycle = 0;
    for (int f = 0; f < fieldNum; f++) stepsPerCycle += mcField[f]->length();

    rando = rando0;
  }

  ~MetroMonteCarlo() {
    // Note that the field objects themselves are owned by DiffusionFusion.
    // We just deallocate our arrays of pointers to them.
    delete[] mcField;
    delete[] mcFieldDesc;
    // rando is owned by the caller.
  }

  int getLocalFieldNum() const { return mcLocalFields.length(); }
  int getFieldNum() const { return fieldNum; }
  IndexList getFieldIds() const {
    IndexList ret;
    for (int f = 0; f < fieldNum; f++) ret.add(mcFieldDesc[f]->id);
    return ret;
  }

  int getStepsPerCycle() const {
    return stepsPerCycle;
  }
  void setStepsPerCycle(int spc) {
    stepsPerCycle = spc;
  }

  // Pick a field node at random and twiddle it.
  TrialMove trialMove(bool local) {
    bool looking = true;
    double jump;
    double currVal;
    trialField = 0;

    while (looking) {
      // Pick the field randomly.
      if (fieldNum != 1) {
	if (local) trialField = mcLocalFields.get( int(rando->uniform()*mcLocalFields.length()) );
	else trialField = int(rando->uniform()*fieldNum);
      }
      const FieldDesc* fd = mcFieldDesc[trialField];
      
      // Pick the field node randomly.
      if (fd->fixedFile.length() == 0)
	// All nodes are free.
	trialNode = int(rando->uniform()*mcField[trialField]->length());
      else {
	// Free nodes are specified in freeNodeList.
	int trialInd = int(rando->uniform()*fd->freeNodeList.length());
	trialNode = fd->freeNodeList.get(trialInd);
      }
      
      // Choose the move.
      jump = fd->step*rando->student(1.0);

      // The move.
      currVal = mcField[trialField]->get(trialNode);
      trialVal = currVal + jump;
      
      // We immediately discard moves that violate the field limits.
      if (fd->hasMin && trialVal < fd->minVal) continue;
      if (fd->hasMax && trialVal > fd->maxVal) continue;
      looking = false;
    }

    // Make the trial move.
    trialLastVal = currVal;
    mcField[trialField]->set(trialNode, trialVal);

    Vector3 trialPos(0.0);
    // Limitation: Only the first three dimensions of the position are kept.
    trialPos.x = mcField[trialField]->nodePos(trialNode, 0);
    if ( mcField[trialField]->dimensions() >= 2 )  
      trialPos.y = mcField[trialField]->nodePos(trialNode, 1);
    if ( mcField[trialField]->dimensions() >= 3 )  
      trialPos.z = mcField[trialField]->nodePos(trialNode, 2);

    return TrialMove(mcFieldDesc[trialField]->id, trialNode, trialVal, trialLastVal, trialPos);
  }

  void reject() {
    // Revert the move.
    mcField[trialField]->set(trialNode, trialLastVal);
  }

  bool metropolis(double deltaCost) {
    if (deltaCost <= 0.0) return true;
    return (rando->uniform() < exp(-deltaCost));
  }

private:
  // Don't permit.
  MetroMonteCarlo();
  MetroMonteCarlo(const MetroMonteCarlo&);
  void operator=(const MetroMonteCarlo&);
};
#endif
