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

// $Header: /private-cvsroot/libraries/arguments/src/arguments.h,v 1.5 2007-02-07 15:12:42 jharlap Exp $

#ifndef ARGUMENTS_H
#define ARGUMENTS_H

#include <string>
#include <map>
#include <vector>

/*! \addtogroup arguments The arguments parser*/
/*@{*/
using namespace std;

/*! A class that handles commandline argument parsing.
 *
 * An Argument is a mandatory argument - ie, a token specified by the
 * user. Options are optional arguments, always prefixed by one or
 * more characters.  For example, a common Option to add is '-help'.
 * Options can have Arguments of their own, such as '-size 3'.
 */

class Arguments {
  string optionPrefix;
  string programName;
  string programDescription;
  
 public:
  class Option;

  /*! A class that holds a single mandatory argument
   *
   * A mandatory argument on the command line, or a mandatory argument
   * to an Option.
   */
  class Argument {
    friend class Arguments;
    friend class Option;
    
    string argName;
    string argDescription;
    string argValue;
    string argDef;
    
  public:
    /*! Argument constructor.
     *
     * sets the name, description and default value of an argument.
     *      
     * \param name the name of the argument
     * \param description a brief description
     * \param def the default value if none is specified
     */
    Argument(string name, string description="", string def="");
  };
  
  typedef vector<Argument> ArgumentVector;
  
  /*! A class that holds a single option.
   *
   * Each option is one optional argument on the command line, and may
   * contain 0 or more mandatory Arguments.
   */
  class Option {
    friend class Arguments;
    
    string optName;
    string optDescription;
    ArgumentVector optArguments;
    bool optIsSet;

  public:
    /*! Option constructor.
     *
     * simply sets the name and description of the new Option.
     *
     * \param name name of the option (for the option "-help", use simply "help" here)
     * \param description a brief description of the option
     */
    Option(string name, string description);

    /*! Adds an Argument to the Option.
     *
     * Forces the Option to be followed by an argument, such as in
     * "-size 3" (where '3' is the Argument).
     *
     * \param name name of the argument
     * \param description a brief description of the argument
     * \param def the default value
     */
    void addArgument(string name, string description="", string def="");

    /*! [] operator for accessing arguments of the Option
     *
     * Provides an associative container type of interface to
     * Arguments of an Option.
     *
     * \param name the name of the Argument to be accessed
     *
     * \return the value of the Argument requested
    */
    string &operator[](const char * name);

    /*! bool operator determines if an option has been set.
     *
     * \return true if the Option was passed on the command line, false otherwise
    */
    operator bool();

    /*! Gets the name of the Option.
     *
     * \return the name of the option
     */
    string getName();
  };

 public:
  /*! A type to hold the leftover arguments from the command line */
  typedef vector<string> LeftoversVector;

 private:
  typedef map<string, Option> OptionMap;

  OptionMap argsOptions;
  ArgumentVector argsArguments;
  static Option unknownOption;
  static Argument unknownArgument;
  LeftoversVector argsLeftovers;


 public:
  /*! Prints out a help/usage screen.
   *
   * Generates a help/usage screen from the children options and
   * arguments of the Arguments object.
   */
  void usage();

  /*! Adds an option.
   *
   * \param name the name of the option
   * \param description the description of the option
   */
  void addOption(string name, string description="");

  /*! Adds an option.
   *
   * \param option the Option object to be added.
   */
  void addOption(Option &option);

  /*! Adds an argument.
   * 
   * \param name the name of the argument
   * \param description a brief description of the argument
   * \param def the default value of the argument
   */
  void addArgument(string name, string description="", string def="");

  /*! Parses the command line
   *
   * Parses the command line (simply pass forward argc and argv) using
   * the existing tree of Options and Arguments.
   *
   * \param argc the number of arguments on the command line
   * \param argv an array of strings, with length equal to argc
   *  
   * \return true if parse() succeeded, false if something went wrong
   */
  bool parse(int argc,char *argv[]);

  /*! [] operator for accessing mandatory arguments.
   *
   * Provides an associative container type of interface to
   * Arguments.
   *
   * \param name the name of the Argument to be accessed
   *
   * \return the value of the Argument requested
   */
  string &operator[](const char* name);

  /*! Gets an Option object by name.
   *
   * \param name the name of the Option wanted
   *
   * \return a reference to the Option
   */
  Option &getOption(const char* name);

  /*! Gets the leftover arguments.
   *
   * The vector of leftover arguments is essentially a vector<string>,
   * and is generally useful if a program can operate a variable
   * number of files.
   *
   * \return the leftover arguments
   */
  LeftoversVector &getLeftovers();

  /*! Arguments constructor.
   *
   * Does nothing except initialize the program name, description and option prefix.
   *
   * \param name the name of the program
   * \param description a brief description of the program
   * \param prefix the prefix to use for Options, defaults to "-"
   */
  Arguments(string name, string description="", string prefix="-");
};
/*@}*/

#endif // ARGUMENTS_H
