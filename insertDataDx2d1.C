///////////////////////////////////////////////////////////////////////
// Author: Jeff Comer <jeffcomer at gmail>

#include "useful.H"
#include "PiecewiseBicubic.H"

int countValueLines(const char* fileName) {
  FILE* inp = fopen(fileName, "r");
  char line[STRLEN];
  int count = 0;

  while (fgets(line, STRLEN, inp) != NULL) {
    // Ignore comments.
    int len = strlen(line);
    if (line[0] == '#') continue;
    if (len < 2) continue;
    count++;
  }
  fclose(inp);

  return count;
}

int readValueLines(const char* fileName, double* val) {
  FILE* inp = fopen(fileName, "r");
  char line[STRLEN];
  int count = 0;

  while (fgets(line, STRLEN, inp) != NULL) {
    // Ignore comments.
    int len = strlen(line);
    if (line[0] == '#') continue;
    if (len < 2) continue;

    val[count] = strtod(line, NULL);
    count++;
  }
  fclose(inp);

  return count;
}

////////////////////////////////////////////////////////////////////////////////
// Driver
int main(int argc, const char* argv[]) {
  if (argc != 4) {
    printf("Usage: %s dataFile inDxFile outDxFile\n", argv[0]);
    exit(0);
  }

  const char* dataFile = argv[1];  
  const char* dxFile = argv[2];
  const char* outFile = argv[argc-1];

  PiecewiseBicubic field(dxFile, false);
  int fieldN = field.length();
  int dataN = countValueLines(dataFile);
  
  if (fieldN != dataN) {
    fprintf(stderr, "ERROR dx `%s' and file `%s' do not have the same number of elements (%d and %d).\n", dxFile, dataFile, fieldN, dataN);
    exit(-1);
  }
  
  double* val = new double[dataN];
  readValueLines(dataFile, val);
  for (int i = 0; i < dataN; i++) field.set(i,val[i]);
  
  field.write(outFile);
  
  return 0;
}

