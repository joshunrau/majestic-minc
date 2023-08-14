// some simple statistics of a vertstats file

#include "mniVertstatsFile.h"
#include "mniVertstatsMath.h"
#include <iostream>
#include <vector>

extern "C" {
#include <ParseArgv.h>
}

// argument defaults
char *columnName = NULL;

// argument handling
ArgvInfo argTable[] = 
  { {NULL, ARGV_HELP, (char *) NULL, (char *) NULL,
     "File input options" },
    {"-column", ARGV_STRING, (char *) 1, (char *) &columnName,
     "Name of the column to load. [Default: load first column]"},
    {NULL, ARGV_END, NULL, NULL, NULL }
  };

using namespace std;

int main( int argc, char *argv[] ) {
  string infile;
  
  string usage = "vertstats_stats [options] filename.vertstats";

  if (ParseArgv( &argc, argv, argTable, 0 ) || (argc != 2) ) {
    cerr << "USAGE: " << usage << endl;
    exit(1);
  }

  infile = argv[1];

  mniVertstatsFile f(infile);

  if (columnName == NULL) {
    // output stats for all columns
    for (int i=0; i < f.getNumColumns(); i++) {
      std::vector<string> header = f.getDataHeader();
      mniVectorStats s(f.getDataColumn(i));
      cout << header[i] << ": " << endl << s;
    }
  }
  else {
    // output stats for one column only
    mniVectorStats s(f.getDataColumn(columnName));
    cout << columnName << ": " << endl << s;
  }

  return(0);
}

  

