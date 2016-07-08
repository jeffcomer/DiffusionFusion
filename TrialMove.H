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
// A struct with all public members for holding a descriptor of a node.
// Author: Jeff Comer <jeffcomer at gmail>
#include "useful.H"

#ifndef TRIALMOVE_H
#define TRIALMOVE_H

struct TrialMove {
  int fieldId;
  int node;
  double trialVal, lastVal;
  Vector3 pos;
  TrialMove() {}
  TrialMove(int id, int n, double tv, double lv, Vector3 pos0) :
    fieldId(id), node(n), trialVal(tv), lastVal(lv), pos(pos0) {
   }
};
#endif
