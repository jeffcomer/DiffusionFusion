///////////////////////////////////////////////////////////////////////
// Author: Jeff Comer <jeffcomer at gmail>

#include "useful.H"
#include "PiecewiseBicubic.H"

int main(int argc, char* argv[]) {
  if (argc != 3) {
    printf("Usage: %s inDxFile outFile\n", argv[0]);
    exit(0);
  }
  const char* inFile = argv[1];
  const char* outFile = argv[argc-1];

  PiecewiseBicubic field(inFile, false);
  const double x0 = field.getX0();
  const double dx = field.getDx();
  const int nx = field.getNx();
  const double y0 = field.getY0();
  const double dy = field.getDy();
  const int ny = field.getNy();
  
  // Open the file.
  FILE* out = fopen(outFile,"w");
  if (out == NULL) {
    printf("ERROR Couldn't open file %s.\n",outFile);
    return false;
  }
  
  // Write the header.
  fprintf(out, "# 2\n");
  fprintf(out, "# %.15g %.15g %d 0\n", x0, dx, nx);
  fprintf(out, "# %.15g %.15g %d 0\n\n", y0, dy, ny);
  
  // Write the data.
  for (int ix = 0; ix < nx; ix++) {
    for (int iy = 0; iy < ny; iy++) {
      double val = field.get(field.index(ix,iy));
      fprintf(out," %.15g %.15g  %.15g\n", field.getPosX(ix), field.getPosY(iy), val);
    }
    // An extra newline between different x values.
    fprintf(out,"\n");
  }
  fclose(out);

  return 0;
}
