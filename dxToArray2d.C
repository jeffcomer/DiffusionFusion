///////////////////////////////////////////////////////////////////////
// Author: Jeff Comer <jeffcomer at gmail>

#include "useful.H"
#include "PiecewiseBicubic.H"

int main(int argc, char* argv[]) {
  if (argc != 3) {
    printf("Usage: %s inDxFile outPrefix\n", argv[0]);
    exit(0);
  }
  
  PiecewiseBicubic field(argv[1], false);
  field.writeArray(argv[argc-1]);
  
  return 0;
}
