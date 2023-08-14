// Extracts a single column from a vertex file

#include "mniVertstatsFile.h"
#include <iostream>
#include <vector>

int main(int argc, char *argv[]) {
  if (argc != 3 ) {
    cerr << "USAGE: vertstats_extract filename column_name" << endl;
    return(1);
  }

  mniVertstatsFile f(argv[1]);

  vertexColumn column = f.getDataColumn(argv[2]);
  for (int i=0; i < f.getNumRows(); i++) {
    cout << column[i] << endl;
  }
  return(0);
}
