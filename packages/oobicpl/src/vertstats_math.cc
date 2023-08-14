#include "mniVertstatsFile.h"
#include "mniVertstatsMath.h"

extern "C" {
#include <ParseArgv.h>
// JPL: values.h does not exist on OS X, so here I'm adding a hack and
// defining the necessary values myself
//#include <values.h>
#include <limits.h>
#include <float.h>
#define MAXSHORT   SHRT_MAX
#define MAXINT     INT_MAX 
#define MAXDOUBLE  DBL_MAX

#include <time_stamp.h>
}

#include <iostream>


using namespace std;

#define DEFAULT_DBL DBL_MAX
#define TRUE 1
#define FALSE 0

typedef enum {
  ADD_OP, SUB_OP, MULT_OP, DIV_OP, SEG_OP, ABS_OP, NORM_OP, UNSPECIFIED_OP
} Operation;

// argument variables
int copyHeader = FALSE;
int oldStyle = FALSE;

double constant = DEFAULT_DBL;
double constant2[2] = {DEFAULT_DBL, DEFAULT_DBL};
double maskRange[2] = {DEFAULT_DBL, DEFAULT_DBL};

char *columnName = NULL;

Operation operation = UNSPECIFIED_OP;

// argument table

ArgvInfo argTable[] = 
  { {NULL, ARGV_HELP, (char *) NULL, (char *) NULL,
     "File output options" },
    //    {"-copy_header", ARGV_CONSTANT, (char *) TRUE, (char *) &copyHeader,
    //     "Copy the header from the first input file" },
    {"-old_style_file", ARGV_CONSTANT, (char *) TRUE, (char *) &oldStyle,
     "Write out using the old format" },
    {NULL, ARGV_HELP, (char *) NULL, (char *) NULL,
     "File input options" },
    {"-column", ARGV_STRING, (char *) 1, (char *) &columnName,
     "Name of the column to load of input files [Default: load first column]"},
    {"-mask_range", ARGV_FLOAT, (char *) 2, (char *) &maskRange,
     "Only apply operation to values in this range." },
    {NULL, ARGV_HELP, (char *) NULL, (char *) NULL, 
     "Specification of Constants" },
    {"-const", ARGV_FLOAT, (char *) 1, (char *) &constant,
     "Specify a constant argument" },
    {"-const2", ARGV_FLOAT, (char *) 2, (char *) &constant2,
     "Specify two constant arguments" },
    {NULL, ARGV_HELP, (char *) NULL, (char *) NULL,
     "Mathematical operations" },
    {"-mult", ARGV_CONSTANT, (char *) MULT_OP, (char *) &operation,
     "Multiply two vectors or a vector and a constant" },
    {"-add", ARGV_CONSTANT, (char *) ADD_OP, (char *) &operation, 
     "Add two vectors or a vector and a constant" },
    {"-div", ARGV_CONSTANT, (char *) DIV_OP, (char *) &operation,
     "Divide two vectors or a vector and a constant" },
    {"-sub", ARGV_CONSTANT, (char *) SUB_OP, (char *) &operation,
     "Subtract two vectors or a vector and a constant" },
    {"-seg", ARGV_CONSTANT, (char *) SEG_OP, (char *) &operation,
     "Segment a vector using range specified by -const2. Segmented areas take on value of 1 unless specified differently using -const" },
    {"-abs", ARGV_CONSTANT, (char *) ABS_OP, (char *) &operation,
     "Take absolute value of a vector component-wise" },
    {"-norm", ARGV_CONSTANT, (char *) NORM_OP, (char *) &operation,
     "Normalises a vector by dividing each element by the mean" },

    {NULL, ARGV_END, NULL, NULL, NULL}
  };
     
// Main programme
int main(int argc, char *argv[]) {

  std::vector<string> infiles;
  string outfile;
  string argString;
  unsigned int nfiles;
  std::vector<vertexColumn> input;
  vertexColumn result;
  bool oldStyleBool;

  // to be used for the history string
  argString = time_stamp(argc, argv);

  if (ParseArgv( &argc, argv, argTable, 0) || (argc < 2)) {
    cerr << "ERROR: bad arguments" << endl;
    exit(1);
  }

  // number of files specified
  nfiles = argc - 2;

  // outfile is always last
  outfile = argv[argc-1];

  // load the vertex columns
  for (int i=1; i <= nfiles; i++) {
    mniVertstatsFile tmp(argv[i]);
    if (columnName != NULL) {
      try {
        input.push_back(tmp.getDataColumn(columnName));
      }
      catch(const mniVertstatsFile::InvalidColumnError e) {
        cerr << "ERROR: column name " << columnName << " does not "
             << "exist in file " << argv[i] << endl;
        exit(1);
      }
    }
    else { // no column name specified
      input.push_back(tmp.getDataColumn(0));
    }
  }

  // determine the type of file to be written out
  if (oldStyle == TRUE) {
    oldStyleBool = true; 
  }
  else {
    oldStyleBool = false;
  }

  // perform the desired operation
  switch (operation) {
  case ADD_OP: // addition
    if (constant == DEFAULT_DBL && nfiles != 2) {
      cerr << "ERROR: for -add you must specify either a constant "
           << "or two input files." << endl;
      exit(1);
    }
    else if (constant == DEFAULT_DBL) {
      result = vectorAdd(input[0], input[1]);
    }
    else {
      result = vectorAdd(input[0], constant);
    }
    break;
  
  case MULT_OP: // multiplication
    if (constant == DEFAULT_DBL && nfiles !=2) {
      cerr << "ERROR: for -mult you must specify either a constant "
           << "or two input files." << endl;
      exit(1);
    }
    else if (constant == DEFAULT_DBL) {
      result = vectorMult(input[0], input[1]);
    }
    else {
      result = vectorMult(input[0], constant);
    }
    break;

  case DIV_OP: // division
    if (constant == DEFAULT_DBL && nfiles != 2) {
      cerr << "ERROR: for -div option you must specify either a constant "
           << "or two input files." << endl;
      exit(1);
    }
    else if (constant == DEFAULT_DBL) {
      result = vectorDiv(input[0], input[1]);
    }
    else {
      result = vectorDiv(input[0], constant);
    }
    break;

  case SUB_OP: // subtraction
    if (constant == DEFAULT_DBL && nfiles != 2) {
      cerr << "ERROR: for -sub option you must specify either a constant "
           << "or two input files." << endl;
      exit(1);
    }
    else if (constant == DEFAULT_DBL) {
      result = vectorSub(input[0], input[1]);
    }
    else {
      result = vectorSub(input[0], constant);
    }
    break;

 case SEG_OP: // segmentation
   if (constant2[0] == DEFAULT_DBL || constant2[1] == DEFAULT_DBL) {
     cerr << "ERROR: for -seg option you must specify a range using "
          << "-const2." << endl;
   }
   else {
     if (constant == DEFAULT_DBL) {
       // segmented areas have value of 1
       result = vectorSeg(input[0], constant2[0], constant2[1]);
     }
     else {
       // output takes on value other than one for segmented areas
       result = vectorSeg(input[0], constant2[0], constant2[1], constant);
     }
   }
   break;

  case ABS_OP: // absolute value
    result = vectorAbsolute(input[0]);
    break;

  case NORM_OP: // normalise
    result = vectorNormalise(input[0]);
    break;
  }

  // write the output file
  mniVertstatsFile output;
  output.putDataColumn(result, "Column0");
  output.addToHistory(argString);
  output.writeFile(outfile, oldStyleBool);
  return 0;
}

