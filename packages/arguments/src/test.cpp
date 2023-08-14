#include <arguments.h>
#include <iostream>

int main(int argc, char *argv[]) {
  Arguments cArg("test", "program description");

  cArg.addArgument("file1", "file1 description");
  cArg.addOption("opt", "first option");
  cArg.addOption("short", "make shorts");
  cArg.addOption("byte", "make bytes");
  cArg.addArgument("outputFile", "file to create");

  cArg.addOption("help", "show this help");

  Arguments::Option cOpt("inclusion", "defines polygon inclusion");
  cOpt.addArgument("type", "all or some", "all");
  cArg.addOption(cOpt);

  if(!cArg.parse(argc, argv)) {
    return 1;
  }
  
  if(cArg.getOption("help")) {
    cArg.usage();
    return 0;
  }
  
  if(cArg.getOption("short")) cout << "short is set" << endl;
  else cout << "short is not set" << endl;

  if(cArg.getOption("inclusion")) cout << "inclusion is set to " << cArg.getOption("inclusion")["type"] << endl;

  cout << "outputFile: " << cArg["outputFile"] << endl;

  cout << "unknown argument: " << cArg["notDefinedAbcd"] << endl;
  cout << "unknown option: " << cArg.getOption("notDefined").getName() << endl;

  Arguments::LeftoversVector leftovers = cArg.getLeftovers();

  if(leftovers.size() > 0) {
    cout << "Leftovers: " << endl;
    for(Arguments::LeftoversVector::iterator it = leftovers.begin(); it != leftovers.end(); ++it) {
      cout << *it << endl;
    }
  }
  
  return 0;
}
