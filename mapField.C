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
