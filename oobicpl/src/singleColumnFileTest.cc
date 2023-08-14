// test whether a single column file can be correctly loaded

#include "mniVertstatsFile.h"
#include <iostream>

int main(int argc, char *argv[]) {

  mniVertstatsFile f("single_column_file.txt");
  
  if (f.getNumColumns() != 1) {
    cerr << "Returned incorrect number of columns: should have been 1, was "
         << f.getNumColumns() << endl;
    return(1);
  }

  vertexColumn test = f.getDataColumn(0);
  if (test[1] > 2.82 && test[1] < 2.80) {
    cerr << "Value returned incorrect. Should have been 2.814, was "
         << test[1] << endl;
    return(1);
  }
  return(0);
}
