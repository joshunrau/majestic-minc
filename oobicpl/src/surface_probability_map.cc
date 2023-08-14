#include <mniVertstatsFile.h>
#include <mniVertstatsMath.h>

#include <iostream>
#include <map>
#include <string>
#include <vector>

// ParseArgv inclusion and related defines
extern "C" {
#include <ParseArgv.h>
#include <time_stamp.h>
}

#define TRUE 1
#define FALSE 0

// argument defaults
int oldStyle = FALSE;
char *columnName = NULL;

ArgvInfo arbTable[] = 
  { {NULL, ARGV_HELP, (char *) NULL, (char *) NULL,
     "File output options" },
    //    {"-old_style_file", ARGV_CONSTANT, (char *) TRUE, (char *) &oldStyle,
    //    "Write out using old format" },
    {NULL, ARGV_HELP, (char *) NULL, (char *) NULL,
     "File input options" },
    {"-column", ARGV_STRING, (char *) 1, (char *) &columnName,
     "Name of the column to load for all input files [Default: load first column" },
    {NULL, ARGV_END, NULL, NULL, NULL}
  };

using namespace std;

int main (int argc, char **argv) {
  std::vector<vertexColumn> input;
  string               outfile;
  string               argString;
  unsigned int         nFiles;
  
  // create the output file
  mniVertstatsFile out;
  out.addToHistory(argc, argv);

  if (ParseArgv( &argc, argv, arbTable, 0 ) || (argc < 3)) {
    cerr << "USAGE: " << argv[0] 
	 << " [options] <infile_1> <infile_n> <outfile>" << endl;
    return 0;
  }

  // number of files specified
  nFiles = argc - 2;

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

  
  vertexColumn commonLabel( input[0] );
  vertexColumn probability( input[0] );

  // the actual work begins here
  for (int v1=0; v1 < commonLabel.size(); ++v1) {
    map<unsigned int, unsigned int> labels;
    // create an entry for each label present and record number of times
    // that label is used at this vertex
    for (int v2=0; v2 < input.size(); ++v2) {
      labels[input[v2][v1]]++;
    }
    map<unsigned int, unsigned int>::iterator highest = labels.begin();
    map<unsigned int, unsigned int>::iterator cur = labels.begin();
    // check which label has the most common representation
    while (cur != labels.end()) {
      if ((*cur).second > (*highest).second)
	highest = cur;
      cur++;
    }
    // output most common label as well as the percentage of 
    // times that label was used at this vertex
    commonLabel[v1] = (*highest).first;
    probability[v1] = (float)(*highest).second / input.size();
  }

  // output the two columns into the outfile
  out.putDataColumn(commonLabel, "common.label");
  out.putDataColumn(probability, "prob.of.label");
  out.writeFile(outfile);

  return 0;
}

