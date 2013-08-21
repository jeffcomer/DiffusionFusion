///////////////////////////////////////////////////////////////////////
// Author: Jeff Comer <jeffcomer at gmail>

#include "useful.H"
#include "PiecewiseCubic.H"

int main(int argc, char* argv[]) {
  if (argc != 3) {
    printf("Usage: %s inFile outPrefix\n", argv[0]);
    exit(0);
  }
  
  PiecewiseCubic field(argv[1], false);
  field.writeIntegral(argv[argc-1], -1.0);
  
  return 0;
}
