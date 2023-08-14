#include "mniVertstatsFile.h"
#include "mniVertstatsMath.h"

extern "C" {
#include <ParseArgv.h>
#include <time_stamp.h>
}

#include <iostream>
#include <math.h>

char *columnName = NULL;

ArgvInfo argTable[] = 
  { {NULL, ARGV_HELP, (char *) NULL, (char *) NULL,
     "File input options" },
    {"-column", ARGV_STRING, (char *) 1, (char *) &columnName,
     "Name of the column to load for all input files [Default: load first column" },
    {NULL, ARGV_END, NULL, NULL, NULL}
  };

// Main programme 
int main(int argc, char *argv[]) {
  std::vector<string> infiles;
  string outfile;
  vertexMatrix input;
  int nFiles;

  // create the output file
  mniVertstatsFile out;
  out.addToHistory(argc, argv);

  if (ParseArgv( &argc, argv, argTable, 0 ) || (argc < 3)) {
    cerr << "USAGE: " << argv[0]
	 << " [options] <infile_1> <infile_n> <outfile>" << endl;
    return 0;
  }

  // number of files specified
  nFiles = argc -2;
  
  // outfile is always last
  outfile = argv[argc-1];

  // load the input files
  for (int i=1; i <= nFiles; i++) {
    mniVertstatsFile tmp(argv[i]);
    cout << "Loading: " << argv[i] << endl;
    if (columnName != NULL) {
      try {
	input.push_back(tmp.getDataColumn(columnName));
      }
      catch(const mniVertstatsFile::InvalidColumnError e) {
	cerr << "ERROR: columnn name " << columnName << " does not "
	     << "exist in file " << argv[i] << endl;
	return 0;
      }
    }
    else { // no column name specified
      input.push_back(tmp.getDataColumn(0));
    }
  }

  vertexColumn means( input[0] );
  vertexColumn sds( input[0] );
  vertexColumn cov( input[0] );

  // the actual work begins here
  for (int v1=0; v1 < means.size(); ++v1) {

    // compute mean
    float total = 0;
    for (int v2=0; v2 < nFiles; ++v2) {
      total += input[v2][v1];
    }
    means[v1] = total / nFiles;

    // compute standard deviation
    total = 0;
    for (int v2=0; v2 < nFiles; ++v2) {
      total += pow(input[v2][v1] - means[v1], 2);
    }
    sds[v1] = sqrt(total / (nFiles - 1));

    // compute the coefficient of variation
    cov[v1] = sds[v1] / means[v1];
  }

  out.putDataColumn(means, "mean");
  out.putDataColumn(sds, "standard.deviation");
  out.putDataColumn(cov, "coefficient.of.variation");
  out.writeFile(outfile);

  return 0;
}

