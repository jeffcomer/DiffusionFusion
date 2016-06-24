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
#ifndef EVENT_H
#define EVENT_H

struct Event {
public: 
  static const int varMax = 8;
  double var[varMax]; // variable value
  double del[varMax]; // difference in variable
  double del0[varMax]; // last difference in variable
  int serial; // where does the event lie in the input file?
  int bias; // link to bias field
  int group; // events can be assigned group numbers
  
  // Organizational parameters put in an unset state.
  Event() : serial(-1), bias(-1), group(0) {}
};

#endif
