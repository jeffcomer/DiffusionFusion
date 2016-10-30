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
// Author: Jeff Comer <jeffcomer at gmail>
// A TrajCostComputer descriptor

#ifndef TRAJCOSTDESC_H
#define TRAJCOSTDESC_H

#include "Field.H"
#include "Event.H"

struct TrajCostDesc {
public:
  Field** fieldList; // Pointer to array of fields.
  IndexList fieldSel; // Indices of fields associated with this TrajCostComputer
  Event* event; // Pointer to array of events.
  int eventNum; // Total number of events.
  double kbt; // Thermal energy
  int leastLocal; // Index of the least local field in fieldList
  int group; // TrajCostComputer applies to events with this group (-1 means all events).
  int dimension; // Dimension along which we compute the displacement
  double timestep; // Timestep of Smoluchowski integration
  double maxHop; // Largest number of nodes hopped in a Smoluchowski solution (time delta t)
  double weight; // Multiply the cost by this factor

  TrajCostDesc(Field** fieldList0, Event* event0, int eventNum0) :
    fieldList(fieldList0), event(event0), eventNum(eventNum0), kbt(1.0), 
    leastLocal(-1), group(0), dimension(-1), timestep(-1.0), maxHop(-1.0), weight(1.0)
  {
  }
private:
  // Force use of the constructor above.
  TrajCostDesc() { }
};

#endif
