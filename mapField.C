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

#include "useful.H"
#include "PiecewiseCubic.H"

int main(int argc, char* argv[]) {
  if (argc != 4) {
    printf("Usage: %s nodeFile valFile outFile\n", argv[0]);
    exit(0);
  }
  
  PiecewiseCubic node(argv[1], false);
  PiecewiseCubic val(argv[2], false);

  Field* out = node.map( static_cast<const Field*>(&val) );
  out->write(argv[argc-1]);
  delete out;
  
  return 0;
}
