// to test whether vertstats files can be read correctly

#include "mniVertstatsFile.h"
#include <iostream>

int main(int argc, char *argv[]) {
  
  // open the test file
  //  mniVertstatsFile f("/home/bic/jason/lmu/analysis/all.vertstats");
  mniVertstatsFile f("new_style_file.vertstats");
  try { 
    vertexColumn test3 = f.getDataColumn("slope.1");
  }
  catch (mniVertstatsFile::InvalidColumnError) {
    cerr << "invalid row" << endl;
    return(1);
  }
  return(0);
}

