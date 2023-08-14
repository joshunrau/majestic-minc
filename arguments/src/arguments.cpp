/*
  libarguments provides simple command line argument parsing for C++ programs
  Copyright (C) 2004  Jonathan Harlap
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include "arguments.h"
#include <iostream>

#ifndef lint
static char rcsid[] = "$Header: /private-cvsroot/libraries/arguments/src/arguments.cpp,v 1.4 2004-09-05 13:52:57 jharlap Exp $";
#endif

using namespace std;

/**************************************
 ** Arguments class                  **
 **************************************/
Arguments::Option Arguments::unknownOption("UNKNOWN", "UNKNOWN");
Arguments::Argument Arguments::unknownArgument("UNKNOWN", "UNKNOWN", "UNKNOWN");

/*! Arguments constructor
 *
 * does nothing other than initialize the program name, description and option prefix
 *  
 * \param name the name of the program
 * \param description a brief description of the program, defaults to ""
 * \param prefix the prefix to use for options, defaults to "-"
 */
Arguments::Arguments(string name, string description, string prefix)
  : programName(name),
    programDescription(description),
    optionPrefix(prefix)
{
}

bool Arguments::parse(int argc, char *argv[]) {
  // parse the arguments
  int prefixLength = optionPrefix.length();

  vector<string> unusedArguments;

  // loop over command line tokens - starting at 1 to skip the program name
  for(int tokenIndex = 1; tokenIndex < argc; ++tokenIndex) {
    string token(argv[tokenIndex]);

    // if it's an option
    if(token.substr(0, prefixLength) == optionPrefix) {

      // find the option
      OptionMap::iterator curOpt = argsOptions.find(token.substr(prefixLength));
      if(curOpt == argsOptions.end()) {
        cerr << programName << " error: Unknown option " << token.substr(prefixLength) << "." << endl;
        usage();
        return false;
      }

      curOpt->second.optIsSet = true;

      // if it has arguments
      for(ArgumentVector::iterator curOptArg = curOpt->second.optArguments.begin();
          curOptArg != curOpt->second.optArguments.end();
          ++curOptArg) {
        
        // go to the next token on the command line
        ++tokenIndex;
        
        // make sure we still have a token to access
        if(tokenIndex >= argc) {
          cerr << programName << " error: Too few arguments for option " << token.substr(prefixLength) << "." << endl;
          usage();
          return false;
        }
        
        // set the argument value to the token
        curOptArg->argValue = argv[tokenIndex];
      }
    } else {
      // not an option - must be an argument
      unusedArguments.push_back(token);
    }
  }
  
  // make sure we have enough arguments
  if(unusedArguments.size() < argsArguments.size()) {
    cerr << programName << " error: Too few arguments." << endl;
    usage();
    return false;
  }
  
  // now handle arguments
  ArgumentVector::iterator argsIterator = argsArguments.begin();
  for(vector<string>::iterator unusedArgsIterator = unusedArguments.begin();
      unusedArgsIterator != unusedArguments.end();
      ++unusedArgsIterator) {
    
    if(argsIterator == argsArguments.end()) {
      argsLeftovers.push_back(*unusedArgsIterator);
    } else {
      argsIterator->argValue = *unusedArgsIterator;
      ++argsIterator;
    }
  }
  
  return true;
}

void Arguments::addOption(string name, string description) {
  // add a new option
  argsOptions.insert(pair<string, Option>(name, Option(name, description)));
}

void Arguments::addOption(Option &option) {
  // add a new option
  argsOptions.insert(pair<string, Option>(option.getName(), option));
}

void Arguments::addArgument(string name, string description, string def) {
  // add a new (mandatory) argument
  argsArguments.push_back(Argument(name, description, def));
}

void Arguments::usage() {
  // print out usage
  cerr << "Usage: " << programName;
  cerr << endl;
  
  cerr << "Options:" << endl;
  for(OptionMap::iterator curOption = argsOptions.begin(); curOption != argsOptions.end(); ++curOption) {
    cerr << "\t" << optionPrefix << curOption->first ;
    cerr << "\t" << curOption->second.optDescription;
    cerr << endl;

    // option arguments
    for(ArgumentVector::iterator curArg = curOption->second.optArguments.begin();
        curArg != curOption->second.optArguments.end();
        ++curArg) {
      cerr << "\t  " << curArg->argName << "\t= " << curArg->argDescription;
      if(curArg->argDef != "") cerr << " [Default: " << curArg->argDef << "]";
      if(curArg->argValue != "") cerr << " [Currently: " << curArg->argValue << "]";
      cerr << endl;
    }
  }

  cerr << "Arguments:" << endl;
  for(ArgumentVector::iterator curArg = argsArguments.begin(); curArg != argsArguments.end(); ++curArg) {
    cerr << "\t" << curArg->argName;
    cerr << "\t" << curArg->argDescription;
    if(curArg->argValue != "") cerr << " [Currently: " << curArg->argValue << "]";
    cerr << endl;

  }

  cerr << endl;
  cerr << programDescription << endl;
}

string &Arguments::operator[](const char * name) {
  // fetch a specific argument
  for(ArgumentVector::iterator argIt = argsArguments.begin(); argIt != argsArguments.end(); ++argIt) {
    if(argIt->argName == name)
      return argIt->argValue;
  }

  return unknownArgument.argValue;
}

Arguments::Option &Arguments::getOption(const char * name) {
  // fetch a specific option
  OptionMap::iterator it = argsOptions.find(name);

  if(it != argsOptions.end())
    return it->second;

  return unknownOption;

}

Arguments::LeftoversVector &Arguments::getLeftovers() {
  return argsLeftovers;
}

/**************************************
 ** Arguments::Option class          **
 **************************************/

Arguments::Option::Option(string name, string description)
  : optName(name), optDescription(description), optIsSet(false)
{
}

void Arguments::Option::addArgument(string name, string description, string def) {
  // add a new mandatory argument
  optArguments.push_back(Argument(name, description, def));
}

string& Arguments::Option::operator[](const char * name) {
  // fetch a specific argument from the option
  for(ArgumentVector::iterator argIt = optArguments.begin(); argIt != optArguments.end(); ++argIt) 
  {
    if(argIt->argName == name)
      return argIt->argValue;
  }

  return unknownArgument.argValue;
}

Arguments::Option::operator bool() {
  return optIsSet;
}

string Arguments::Option::getName() {
  return optName;
}

/**************************************
 ** Arguments::Argument class        **
 **************************************/

Arguments::Argument::Argument(string name, string description, string def)
  : argName(name), argDescription(description), argDef(def)
{
}

