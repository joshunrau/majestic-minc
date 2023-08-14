// testing an old style multiple columns file

#include "mniVertstatsFile.h"
#include <iostream>

int main(int argc, char *argv[]) {

  mniVertstatsFile f("multiple_columns_file.vertstats");

  if (f.getNumColumns() != 4) {
    cerr << "Wrong number of columns. Should have been 4, was "
         << f.getNumColumns() << endl;
    return(1);
  }

  vertexColumn test = f.getDataColumn(3);
  if (test[1] > -0.84 && test[1] < -0.86) {
    cerr << "Wrong value returned. Should have been -0.85141, was "
         << test[1] << endl;
    return(1);
  }

  // test to see if retrieval of named columns works:
  vertexColumn test2;
  try {
    test2 = f.getDataColumn("Column3");
  }
  catch (mniVertstatsFile::InvalidColumnError) {
    cerr << "Invalid column: Column3. Column selection does not work." 
         <<  endl;
    return(1);
  }
  if (test2[1] != test[1]) {
    cerr << "Wrong value returned. Should have been -0.85141, was "
         << test2[1] << endl;
    return(1);
  }
  return(0);
}
    
