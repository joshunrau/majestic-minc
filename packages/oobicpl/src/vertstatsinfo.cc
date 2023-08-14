// Provides a bit of info about a vertex file

#include "mniVertstatsFile.h"
#include "mniVertstatsMath.h"
#include <iostream>
#include <vector>
#include <string.h>

string usage = "\n"
"Options:\n"
"  -header [name]\n"
"     Print the contents of the data header 'name'\n"
"  -dataheaders\n"
"     Print the dataheaders\n"
"  -headerstructure\n"
"     Print out the header tree\n"
"  -numrows\n"
"     Print the number of rows in the file\n"
"  -numcolumns\n"
"     Print the number of columns in the file\n"
"  -help\n"
"     Print this help message\n"
"\n"
"USAGE: vertstatsinfo [options] filename.vertstats\n"
"\n";

int main(int argc, char *argv[]) {

  int i;
  std::vector<string> headers;
  std::vector<string>::iterator it;

  // booleans that control what gets printed
  bool dataHeaders = false;
  bool headerStructure = false;
  bool readData = false;
  bool numRows = false;
  bool numColumns = false;
  bool stats = false;
  bool all = true;

  if (argc < 2) {
    cerr << "USAGE: [options] verstatsinfo filename" << endl;
    return(1);
  }

  // argument parsing
  for (i=0; i < argc; i++) {
    if ( strcmp(argv[i], "-help") == 0 ) {
      cerr << usage;
      exit(0);
    }
    else if ( strcmp(argv[i], "-header") == 0 )
      headers.push_back(argv[i+1]);
    else if ( strcmp(argv[i], "-dataheaders") == 0 ) {
      all = false;
      dataHeaders = true;
    }
    else if ( strcmp(argv[i], "-headerstructure") == 0 ) { 
      all = false;
      headerStructure = true;
    }
    else if ( strcmp(argv[i], "-numrows") == 0 ) {
      all = false;
      readData = true;
      numRows = true;
    }
    else if ( strcmp(argv[i], "-stats") == 0 ) {
      all = false;
      readData = true;
      stats = true;
    }
    else if ( strcmp(argv[i], "-numcolumns") == 0 ) {
      all = false;
      numColumns = true;
    }
  }

  if (all) {
    readData = true;
    numColumns = true;
    numRows = true;
    dataHeaders = true;
    headerStructure = true;
  }

  mniVertstatsFile f(argv[argc - 1], readData);
  std::vector<string> header = f.getDataHeader();

  if (headerStructure) {
    cout << "Header structure: " << endl;
    f.printHeaderStructure();
  }

  if (numColumns) 
    cout << "Number of columns: " << f.getNumColumns() << endl;
  if (numRows)
    cout << "Number of rows:    " << f.getNumRows() << endl;

  if (dataHeaders) {
    cout << "Data headers:" << endl;

    for (i=0; i < f.getNumColumns(); i++) {
      cout << "    " << header[i] << endl;
    }
  }

  for (it = headers.begin(); it < headers.end(); it++) {
    cout << "Header: " << *it << ":" << endl;
    cout << f.getHeaderValue(*it) << endl << endl;
  }

  if (stats) {
    for (i=0; i < f.getNumColumns(); i++) {
      mniVectorStats s(f.getDataColumn(i));
      cout << header[i] <<  ": " << endl << s;
    }
    /*
    std::vector<float> test1 = f.getDataColumn(0);
    std::vector<float> test2 = f.getDataColumn(1);
    std::vector<float> test = vectorAdd(test1, test2);
    for (i=0; i < test.size(); i++){
      cout << test1[i] << endl;
      cout << test[i] << endl;
    }
    */
 
  }

  return(0);
}
