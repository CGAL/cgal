#ifndef CGAL_TEST_EXCEPTIONS_H
#define CGAL_TEST_EXCEPTIONS_H

#include <CGAL/Object.h>
#include <CGAL/exceptions.h>
#include <CGAL/Arr_tags.h>

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>

/*
* The test test_traits has a global configuration flag, abort_on_error.
* It determines what happens when an unexpected CGAL assertion, pre-condition,
* post-condition, or error occur. By default abort_on_error is false.
* This means that when an unexpected CGAL assertion, pre-condition, 
* post-condition or error occur, the test does not abort, instead it proceeds
* to the next sub-test.
*
* in general you may test any violation by appending _precondition or 
* _postcondition or _assertion or _warning to the wrappers token in the test 
* input file
*
* the CGAL error and warning handling is set to a failure_handler function 
* that throws a special exceptions, which indicates whether the violation was
* expected or not unexpected. Depending on abort_on_error the right exceptions
* is thrown. the exceptions are caught in perform function.
* so basiclly we have 4 cases:

                         | violation occurred       | violation did 
                         |                          |  not occurred
------------------------------------------------------------------------------
violation is expected    |                          |      fail, if 
(violation appended      |     pass, continue       |  !abort_on_error
to token)                |                          | continue else abort
------------------------------------------------------------------------------
violation is unexpected  | fail, if !abort_on_error |   pass, continue
(use regular tokens)     |   continue else abort    |
------------------------------------------------------------------------------
*/

enum Exception_type {EXPECTED_CONTINUE, EXPECTED_ABORT,
                      UNEXPECTED_CONTINUE, UNEXPECTED_ABORT};
enum Violation_type {NON, PRECONDITION, POSTCONDITION,
                      ASSERTION, WARNING};
enum Read_object_type {POINT, CURVE, XCURVE};

enum Enum_type {NUMBER, SIGN, CURVE_END, BOUNDARY, PARAMETER_SPACE};

std::map<Violation_type,std::string> violation_map;

//Indicates whether the end-of-line has been printed
bool end_of_line_printed;

//indicates to run only one error test or run all the tests
bool abort_on_error;

//indicates if precondition or postcondition or 
//assertion or warning violation is tested
Violation_type violation_tested;

//indicates if precondition or postcondition or 
//assertion or warning is violated
Violation_type violation_occurred;

void print_end_of_line()
{
  std::cout << std::endl;
  end_of_line_printed = true;
}

void print_result(bool result)
{
  if (!result && !end_of_line_printed) print_end_of_line();
    std::cout << ((result) ? "Passed" : "Failed");
  print_end_of_line();
}

std::string global_type;
std::string global_expr;
std::string global_lib;
std::string global_file;
int global_line;
std::string global_msg;

class Test_exception : public CGAL::Failure_exception
{
  private :
   Exception_type e;
  public :
    //constructor
    Test_exception(Exception_type e_tmp) :
      CGAL::Failure_exception(global_lib, global_expr, global_file,
                              global_line, global_msg, global_type)
    {
      e=e_tmp;
    }
    //destructor
    ~Test_exception() throw()
    {}

    Exception_type get()
    {
      return e;
    }
};

//function that throws exceptions with the appropriate data 
//respectivly to the type of the violation

void throw_exceptions(bool q)
{
  print_result(q);
  if (q)//expected violation occurred
    if (abort_on_error)
      throw Test_exception(EXPECTED_ABORT);
    else
      throw Test_exception(EXPECTED_CONTINUE);
  else  //unexpected violation occurred
    if (abort_on_error)
      throw Test_exception(UNEXPECTED_ABORT);
    else
      throw Test_exception(UNEXPECTED_CONTINUE);
}

void did_violation_occur()
{
  if (violation_tested!=NON && violation_tested!=violation_occurred)
  {
    //violation is expected but it did not occur
    throw_exceptions(false);
  }
  //else - taken care by failure_handler function
}

void failure_handler(const char * type, const char * expr, const char * file,
                        int line, const char * msg)
{
  std::string tmp = std::string(file);
  unsigned int loc = tmp.find_last_of("/"); // !!!for windows use back slash !!!to ask efi
  if ((int)loc==-1)
  {
    global_lib=std::string("");
    global_file=tmp;
  }
  else
  {
    global_lib=tmp.substr(0,loc+1);
    global_file=tmp.substr(loc+1);
  }
  global_type = std::string(( !type ? "" : type ));
  global_expr = std::string(( !expr ? "" : expr ));
  global_line = line;
  global_msg = std::string(( !msg ? "" : msg ));
  if (!global_type.compare("precondition"))
  {
    violation_occurred = PRECONDITION;
    throw_exceptions(violation_tested == PRECONDITION);
  }
  else if (!global_type.compare("postcondition"))
  {
    violation_occurred = POSTCONDITION;
    throw_exceptions(violation_tested == POSTCONDITION);
  }
  else if (!global_type.compare("assertion"))
  {
    violation_occurred = ASSERTION;
    throw_exceptions(violation_tested == ASSERTION);
  }
  else if (!global_type.compare("warning"))
  {
    violation_occurred = WARNING;
    throw_exceptions(violation_tested == WARNING);
  }
}

CGAL::Failure_function prev_error_handler;
CGAL::Failure_function prev_warning_handler;

#endif
