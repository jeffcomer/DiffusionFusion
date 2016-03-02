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

int readValueLines(const char* fileName, int column, double* val) {
  FILE* inp = fopen(fileName, "r");
  char line[STRLEN];
  int count = 0;
  int ln = 0;

  while (fgets(line, STRLEN, inp) != NULL) {
    ln++;
    // Ignore comments.
    if (line[0] == '#') continue;

    String s(line);
    // Ignore empty lines.
    s = s.trim();
    if (s.length() == 0) continue;

    // Check for the right number of columns.
    int tokN = s.tokenCount();
    if (tokN <= column) {
      fprintf(stderr, "ERROR line %d of file `%s' has too few columns (%d of %d).\n", ln, fileName, tokN, column);
      exit(-1);
    }

    // Get the desired line.
    String* tokList = new String[tokN];
    s.tokenize(tokList);
    
    val[count] = strtod(tokList[column], NULL);
    count++;
  }
  fclose(inp);

  return count;
}

////////////////////////////////////////////////////////////////////////////////
// Driver
int main(int argc, const char* argv[]) {
  if (argc != 4 && argc != 5) {
    printf("Usage: %s dataFile [dataColumn] inDxFile outDxFile\n", argv[0]);
    exit(0);
  }

  const char* dataFile = argv[1];
  const char* dxFile = argv[argc-2];
  const char* outFile = argv[argc-1];

  // Column from which the data is chosen.
  int dataCol = 0;
  if (argc == 5) dataCol = atoi(argv[2]);

  PiecewiseBicubic field(dxFile, false);
  int fieldN = field.length();
  int dataN = countValueLines(dataFile);
  
  if (fieldN != dataN) {
    fprintf(stderr, "ERROR dx `%s' and file `%s' do not have the same number of elements (%d and %d).\n", dxFile, dataFile, fieldN, dataN);
    exit(-1);
  }
  
  double* val = new double[dataN];
  readValueLines(dataFile, dataCol, val);
  for (int i = 0; i < dataN; i++) field.set(i,val[i]);
  
  field.write(outFile);
  
  return 0;
}

