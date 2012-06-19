#ifndef CGAL_TRAITS_BASE_TEST_H
#define CGAL_TRAITS_BASE_TEST_H

#include <CGAL/basic.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
// #include <cstdlib>

#include <boost/lexical_cast.hpp>

#include <CGAL/exceptions.h>
#include <CGAL/Object.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_enums.h>

#include "IO_test.h"

/* The test test_traits has a global configuration flag, abort_on_error.
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
 * 
 *                          | violation occurred       | violation did 
 *                          |                          |  not occurred
 * ---------------------------------------------------------------------------
 * violation is expected    |                          |      fail, if 
 * (violation appended      |     pass, continue       |  !abort_on_error
 * to token)                |                          | continue else abort
 * ---------------------------------------------------------------------------
 * violation is unexpected  | fail, if !abort_on_error |   pass, continue
 * (use regular tokens)     |   continue else abort    |
 * ---------------------------------------------------------------------------
 */

template <typename T_Traits>
class Traits_base_test : public IO_test<T_Traits> {
protected:
  typedef T_Traits                                      Traits;
  typedef IO_test<Traits>                               IO_test_traits;
  typedef typename IO_test_traits::Point_2              Point_2;
  typedef typename IO_test_traits::X_monotone_curve_2   X_monotone_curve_2;
  typedef typename IO_test_traits::Curve_2              Curve_2;
  
  enum Exception_type {EXPECTED_CONTINUE,
                       EXPECTED_ABORT,
                       UNEXPECTED_CONTINUE,
                       UNEXPECTED_ABORT};
  
  enum Violation_type {NON, PRECONDITION,
                       POSTCONDITION,
                       ASSERTION,
                       WARNING};
  
  enum Read_object_type {POINT, CURVE, XCURVE};

  enum Enum_type {NUMBER, SIGN, CURVE_END, BOUNDARY, PARAMETER_SPACE};

  /*! The input data file of points*/
  std::string m_filename_points;

  /*! The input data file of curves*/
  std::string m_filename_curves;

  /*! The input data file of xcurves*/
  std::string m_filename_xcurves;

  /*! The input data file of commands*/
  std::string m_filename_commands;

  /*! The traits type */
  std::string m_traitstype;

  /*! The container of input points */
  std::vector<Point_2>  m_points;

  /*! The container of input curves */
  std::vector<Curve_2>  m_curves;

  /*! The container of x-monotone curves */
  std::vector<X_monotone_curve_2>  m_xcurves;

  /*! Indicates whether the end-of-line has been printed */
  bool m_eol_printed;

  /*! Indicates whether abort after the first failure or to continue */
  bool m_abort_on_error;

  std::map<Violation_type,std::string> m_violation_map;

  //indicates if precondition or postcondition or 
  //assertion or warning is violated
  Violation_type m_violation_occurred;

  //indicates if precondition or postcondition or 
  //assertion or warning violation is tested
  Violation_type m_violation_tested;
  
  /*! Collect the data depending on obj_t */
  bool read_input(std::ifstream& is, Read_object_type obj_t);

  /*! Execute a command */
  virtual bool exec(std::istringstream& str_stream,
                    const std::string& str_command,
                    bool& result) = 0;
  
  /*! Perform the test */
  bool perform(std::ifstream& is);

  /*! Skip comments */
  void skip_comments(std::ifstream& is, char* one_line);

  std::string remove_blanks(char* str);

  bool get_expected_boolean(std::istringstream& str_stream);

  unsigned int get_expected_enum(std::istringstream& str_stream);

  bool translate_boolean(std::string& str_value);

  unsigned int translate_enumerator(std::string& str_value);

  std::pair<Enum_type, unsigned int>
  translate_int_or_text(std::string& str_value);

  std::pair<Enum_type, unsigned int>
  get_next_input(std::istringstream& str_stream);

  /*! Print curve-end string */
  const char* curve_end_str(CGAL::Arr_curve_end cv_end) const
  {
    return (cv_end == CGAL::ARR_MIN_END) ? "MIN_END" : "MAX_END";
  }
  
  /*! Print an error message */
  void print_error(const std::string& msg)
  {
    std::cerr << "Error: " << msg.c_str() << std::endl;
  }
  
  /*! Print the end-of-line */
  void print_eol()
  {
    std::cout << std::endl;
    m_eol_printed = true;
  }

  /*! Print information */
  void print_info(std::string& info,
                  bool start_line = true, bool end_line = true)
  {
    if (start_line && !m_eol_printed) print_eol();
    std::cout << info.c_str();
    if (end_line) print_eol();
  }

  /*! Print final results */
  void print_result(bool result)
  {
    std::string result_str((result) ? "Passed" : "Failed");
    print_info(result_str, false);
  }

  /*! Print expected answer and real answer */
  void print_answer(const std::string& exp, const std::string& real,
                    const std::string& str)
  {
    print_info(std::string("Expected ").append(str).append(": ").append(exp));
    print_info(std::string("Obtained ").append(str).append(": ").append(real));
  }

  /*! Compare two points */
  bool compare_points(const Point_2& exp_answer, const Point_2& real_answer)
  {
    typename Traits::Equal_2 equal = this->m_traits.equal_2_object();
    if (equal(exp_answer, real_answer)) return true;

    std::string exp_answer_str = boost::lexical_cast<std::string>(exp_answer);
    std::string real_answer_str = boost::lexical_cast<std::string>(real_answer);
    print_answer(exp_answer_str, real_answer_str, "point");
    return false;
  }

  /*! Compare two x-monotone curves */
  bool compare_curves(const X_monotone_curve_2& exp_answer,
                      const X_monotone_curve_2& real_answer)
  {
    typename Traits::Equal_2 equal = this->m_traits.equal_2_object();
    if (equal(exp_answer, real_answer)) return true;

    std::string exp_answer_str = boost::lexical_cast<std::string>(exp_answer);
    std::string real_answer_str = boost::lexical_cast<std::string>(real_answer);
    print_answer(exp_answer_str, real_answer_str, "x-monotone curve");
    return false;
  }

  /*! Compare two unsigned int */
  bool compare(const unsigned int& exp_answer,
               const unsigned int& real_answer,
               const char* str = "result")
  {
    if (exp_answer == real_answer) return true;
    std::string exp_answer_str = boost::lexical_cast<std::string>(exp_answer);
    std::string real_answer_str = boost::lexical_cast<std::string>(real_answer);
    print_answer(exp_answer_str, real_answer_str, str);
    return false;
  }

public:
  /*! Constructor */
  Traits_base_test(int argc, char* argv[]);

  /*! Destructor */
  virtual ~Traits_base_test();

  /*! Entry point */
  bool start();
};

/*!
 * Constructor. 
 * Accepts test data file name.
 */
template <class T_Traits>
Traits_base_test<T_Traits>::Traits_base_test(int argc, char* argv[]) :
  m_eol_printed(true),
  m_abort_on_error(false)       // run all tests
{
  if (argc != 6) {
    print_info(std::string("Usage: ").append(argv[0]).
               append(" points_file xcurves_file curves_file commands_file traits_type_name"));
    return;
  }

  typedef T_Traits Traits;
  m_violation_map[PRECONDITION] = std::string("precondition");
  m_violation_map[POSTCONDITION] = std::string("postcondition");
  m_violation_map[ASSERTION] = std::string("assertion");
  m_violation_map[WARNING] = std::string("warning");
  
  m_filename_points = argv[1];
  m_filename_xcurves = argv[2];
  m_filename_curves = argv[3];
  m_filename_commands = argv[4];
  m_traitstype = argv[5];
}

/*!
 * Destructor. 
 */
template <class T_Traits>
Traits_base_test<T_Traits>::~Traits_base_test()
{
  m_filename_points.clear();
  m_filename_xcurves.clear();
  m_filename_curves.clear();
  m_filename_commands.clear();
  m_points.clear();
  m_curves.clear();
  m_xcurves.clear();
}

/*!
 * Test entry point 
 */
template<class T_Traits>
bool Traits_base_test<T_Traits>::start()
{
  std::ifstream in_pt(m_filename_points.c_str());
  std::ifstream in_xcv(m_filename_xcurves.c_str());
  std::ifstream in_cv(m_filename_curves.c_str());
  std::ifstream in_com(m_filename_commands.c_str());
  if (!in_pt.is_open()) {
    print_error(std::string("cannot open file ").append(m_filename_points));
    return false;
  }
  if (!in_xcv.is_open()) {
    print_error(std::string("cannot open file ").append(m_filename_xcurves));
    return false;
  }
  if (!in_cv.is_open()) {
    print_error(std::string("cannot open file ").append(m_filename_curves));
    return false;
  }
  if (!in_com.is_open()) {
    print_error(std::string("cannot open file ").append(m_filename_commands));
    return false;
  }
  if (!read_input(in_pt, POINT)) {
    in_pt.close(); 
    return false;
  }
  if (!read_input(in_xcv, XCURVE)) {
    in_xcv.close(); 
    return false;
  }
  if (!read_input(in_cv,CURVE)) {
    in_cv.close(); 
    return false;
  }
  if (!perform(in_com)) {
    in_com.close(); 
    return false;
  }
  return true;
}

template <class T_Traits>
bool Traits_base_test<T_Traits>::read_input(std::ifstream& is,
                                       Read_object_type obj_t)
{
  char one_line[1024];
  skip_comments(is, one_line);
  std::istringstream str_stream(one_line, std::istringstream::in);
  try {
    for (int i = 0; !is.eof() ; ++i) {
      switch (obj_t) {
       case POINT :
        m_points.resize(m_points.size()+1);
        if (!read_point(str_stream, m_points[i])) {
          print_error(std::string("failed to read point!"));
          return false;
        }
        break;
       case CURVE :
        m_curves.resize(m_curves.size()+1);
        if (!read_curve(str_stream, m_curves[i])) {
          print_error(std::string("failed to read curve!"));
          return false;
        }
        break;
       case XCURVE :
        m_xcurves.resize(m_xcurves.size()+1);
        if (!read_xcurve(str_stream, m_xcurves[i])) {
          print_error(std::string("failed to read xcurve!"));
          return false;
        }
        break;
      }
      str_stream.clear();
      skip_comments(is, one_line);
      str_stream.str(one_line);
    }
  }
  catch (std::exception e) {
    print_error(std::string("exception!"));
    is.close();
    return false;
  }
  return true;
}

/*!
 * Command dispatcher. Retrieves a line from the input file and performes
 * some action. See comments for suitable function in order to know specific
 * command arguments. 
 */
template <class T_Traits>
bool Traits_base_test<T_Traits>::perform(std::ifstream& is)
{
  bool test_result = true;
  std::cout << "Performing test : traits type is " << m_traitstype
            << ", input files are " << m_filename_points << " "
            << m_filename_xcurves << " " << m_filename_curves << " "
            << m_filename_commands << std::endl;
  m_eol_printed = true;
  char one_line[1024];
  char buff[1024];
  bool abort = false;
  int counter = 0;
  while (!(is.eof() || abort)) {
    skip_comments(is, one_line);
    std::istringstream str_stream(one_line, std::istringstream::in);
    buff[0] = '\0';
    str_stream.getline(buff, 1024, ' ');
    std::string str_command(buff);
    std::size_t location = 0;
    m_violation_occurred = m_violation_tested = NON;

    if ((int)str_command.find("_precondition", 0) != -1) {
      location = str_command.find("_precondition", 0);
      m_violation_tested = PRECONDITION;
    }
    else if ((int)str_command.find("_postcondition",0) != -1) {
      location = str_command.find("_postcondition", 0);
      m_violation_tested = POSTCONDITION;
    }
    else if ((int)str_command.find("_assertion", 0) != -1) {
      location = str_command.find("_assertion", 0);
      m_violation_tested = ASSERTION;
    }
    else if ((int)str_command.find("_warning", 0) != -1) {
      location = str_command.find("_warning", 0);
      m_violation_tested = WARNING;
    }

    counter++;
    std::cout << "case number : " << counter << std::endl;
    if (m_violation_tested != NON) {
#if !defined(CGAL_NDEBUG)
      str_command = str_command.substr(0, location);
      std::cout << "Test " << m_violation_map[m_violation_tested] 
                << " violation : ";
#else
      std::cout << "Skipping condition tests in release mode." << std::endl;
      continue;
#endif
    }
    
    try {
      bool result;
      bool ignore = exec(str_stream, str_command, result);
      if (ignore) continue;
      if ((m_violation_tested != NON) &&
          (m_violation_tested != m_violation_occurred))
      {
        //violation is expected but it did not occur
        result = false;
        if (m_abort_on_error) abort = true;
      }
      print_result(result);
      test_result &= result;
    }
    catch (CGAL::Precondition_exception /* e */) {
      if (m_violation_tested != PRECONDITION) {
        test_result = false;
        if (m_abort_on_error) abort = true;
      }
    }
    catch (CGAL::Postcondition_exception /* e */) {
      if (m_violation_tested != POSTCONDITION) {
        test_result = false;
        if (m_abort_on_error) abort = true;
      }
    }
    catch (CGAL::Warning_exception /* e */) {
      if (m_violation_tested != WARNING) {
        test_result = false;
        if (m_abort_on_error) abort = true;
      }
    }
    catch (CGAL::Assertion_exception /* e */) {
      if (m_violation_tested != ASSERTION) {
        test_result = false;
        if (m_abort_on_error) abort = true;
      }
    }
  }
  return test_result;
}

/*!
 * Skip comments. Comments start with the '#' character and extend to the
 * end of the line
 */
template <class T_Traits>
void Traits_base_test<T_Traits>::skip_comments(std::ifstream& is,
                                               char* one_line)
{
  while (!is.eof()) {
    is.getline(one_line, 1024);
    if (one_line[0] != '#') break;
  }
}

/*!
 */
template <class T_Traits>
std::string Traits_base_test<T_Traits>::remove_blanks(char* str)
{
  std::string result = "";
  bool flag = false;
  //only alphanumeric characters and underscores are allowed
  for (; *str != '\0'; ++str)
  {
    if ((*str >= '0' && *str <= '9') || //digits
        (*str >= 'A' && *str <= 'Z') || //upper case letters
        (*str >= 'a' && *str <= 'z') || //lower case letters
         *str == '_') //underscores
    {
      if (!flag)
        flag=true;
      result += *str;
    }
    if (*str == ' ' && flag)
      break;
  }
  return result;
}

/*!
 */
template <class T_Traits_class>
bool Traits_base_test<T_Traits_class>::
translate_boolean(std::string& str_value)
{
  if (str_value == "TRUE") return true;
  return false;
}

/*!
 */
template <class T_Traits>
unsigned int
Traits_base_test<T_Traits>::translate_enumerator(std::string& str_value)
{
  if (str_value == "LARGER" ) {
    return static_cast<unsigned int>(CGAL::LARGER);
  } else if (str_value == "SMALLER" ) {
    return static_cast<unsigned int>(CGAL::SMALLER);
  } else if (str_value == "EQUAL" ) {
    return static_cast<unsigned int>(CGAL::EQUAL);
  }
  CGAL_error();
  return static_cast<unsigned int>(-220776); // My birthday :-)
}

/*!
 */
template <class T_Traits>
std::pair<typename Traits_base_test<T_Traits>::Enum_type, unsigned int>
Traits_base_test<T_Traits>::translate_int_or_text(std::string& str_value)
{
  if (str_value == "MIN_END" ) {
    return std::pair<enum Enum_type, unsigned int>(CURVE_END,CGAL::ARR_MIN_END);
  } else if (str_value == "MAX_END" ) {
    return std::pair<enum Enum_type, unsigned int>(CURVE_END,CGAL::ARR_MAX_END);
  } else if (str_value == "SMALLER" ) {
    return std::pair<enum Enum_type, unsigned int>
      (SIGN, static_cast<unsigned int>(CGAL::SMALLER));
  } else if (str_value == "EQUAL" ) {
    return std::pair<enum Enum_type,unsigned int>(SIGN,
                                    static_cast<unsigned int>(CGAL::EQUAL));
  } else if (str_value == "LARGER" ) {
    return std::pair<enum Enum_type,unsigned int>(SIGN,
                                    static_cast<unsigned int>(CGAL::LARGER));
  } else if (str_value == "ARR_INTERIOR" ) {
    return std::pair<enum Enum_type,unsigned int>
      (BOUNDARY,static_cast<unsigned int>(CGAL::ARR_INTERIOR));
  } else if (str_value == "BOTTOM_BOUNDARY" ) {
    return std::pair<enum Enum_type,unsigned int>
      (PARAMETER_SPACE,static_cast<unsigned int>(CGAL::ARR_BOTTOM_BOUNDARY));
  } else if (str_value == "TOP_BOUNDARY" ) {
    return std::pair<enum Enum_type,unsigned int>
      (PARAMETER_SPACE,static_cast<unsigned int>(CGAL::ARR_TOP_BOUNDARY));
  } else if (str_value == "LEFT_BOUNDARY" ) {
    return std::pair<enum Enum_type,unsigned int>
      (PARAMETER_SPACE,static_cast<unsigned int>(CGAL::ARR_LEFT_BOUNDARY));
  } else if (str_value == "RIGHT_BOUNDARY" ) {
    return std::pair<enum Enum_type,unsigned int>
      (PARAMETER_SPACE,static_cast<unsigned int>(CGAL::ARR_RIGHT_BOUNDARY));
  } else if (str_value == "INTERIOR" ) {
    return std::pair<enum Enum_type,unsigned int>
      (PARAMETER_SPACE,static_cast<unsigned int>(CGAL::ARR_INTERIOR));
  }
  return std::pair<enum Enum_type,unsigned int>
    (NUMBER,static_cast<unsigned int>(std::atoi(str_value.c_str())));
}

/*!
 */
template <class T_Traits>
bool
Traits_base_test<T_Traits>::
get_expected_boolean(std::istringstream& str_stream)
{
  char buff[1024];
  str_stream.getline( buff, 1024, '.');
  buff[str_stream.gcount()] = '\0';
  std::string str_expres = remove_blanks(buff);
  return translate_boolean(str_expres);
}

/*!
 */
template <class T_Traits>
unsigned int
Traits_base_test<T_Traits>::get_expected_enum(std::istringstream& str_stream)
{
  char buff[1024];
  str_stream.getline(buff, 1024, '.');
  buff[str_stream.gcount()] = '\0';
  std::string str_expres = remove_blanks(buff);
  return translate_enumerator(str_expres);
}

/*!
 */
template <class T_Traits>
std::pair<typename Traits_base_test<T_Traits>::Enum_type, unsigned int>
Traits_base_test<T_Traits>::get_next_input(std::istringstream& str_stream)
{
  char buff[1024];
  do {
    str_stream.getline(buff, 1024, ' ');
  } while (str_stream.gcount() == 1);
  buff[str_stream.gcount()] = '\0';
  std::string str_expres = remove_blanks(buff);
  return translate_int_or_text(str_expres);
}

#endif
