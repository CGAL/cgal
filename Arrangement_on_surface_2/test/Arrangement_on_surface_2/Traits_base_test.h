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
#include <CGAL/use.h>

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

template <typename Geom_traits_T>
class Traits_base_test : public IO_test<Geom_traits_T> {
protected:
  typedef Geom_traits_T                                 Traits;
  typedef IO_test<Traits>                               Base;

  typedef typename Base::Point_2                        Point_2;
  typedef typename Base::X_monotone_curve_2             X_monotone_curve_2;
  typedef typename Base::Curve_2                        Curve_2;

  typedef typename Base::Points_vector                  Points_vector;
  typedef typename Base::Xcurves_vector                 Xcurves_vector;
  typedef typename Base::Curves_vector                  Curves_vector;

  enum Exception_type {EXPECTED_CONTINUE,
                       EXPECTED_ABORT,
                       UNEXPECTED_CONTINUE,
                       UNEXPECTED_ABORT};

  enum Violation_type {NON, PRECONDITION,
                       POSTCONDITION,
                       ASSERTION,
                       WARNING};

  enum Enum_type {NUMBER, SIGN, CURVE_END, BOUNDARY, PARAMETER_SPACE};

  /*! The input data file of commands*/
  std::string m_filename_commands;

  /*! The traits type */
  std::string m_traitstype;

  /*! Indicates whether abort after the first failure or to continue */
  bool m_abort_on_error;

  std::map<Violation_type,std::string> m_violation_map;

  //indicates if precondition or postcondition or
  //assertion or warning is violated
  Violation_type m_violation_occurred;

  //indicates if precondition or postcondition or
  //assertion or warning violation is tested
  Violation_type m_violation_tested;

  /*! Execute a command */
  virtual bool exec(std::istringstream& str_stream,
                    const std::string& str_command,
                    bool& result) = 0;

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
  { return (cv_end == CGAL::ARR_MIN_END) ? "MIN_END" : "MAX_END"; }

  /*! Compare two points */
  bool compare_points(const Point_2& exp_answer, const Point_2& real_answer)
  {
    typename Traits::Equal_2 equal = this->m_geom_traits.equal_2_object();
    if (equal(exp_answer, real_answer)) return true;

    std::string exp_answer_str = boost::lexical_cast<std::string>(exp_answer);
    std::string real_answer_str = boost::lexical_cast<std::string>(real_answer);
    this->print_answer(exp_answer_str, real_answer_str, "point");
    return false;
  }

  /*! Compare two x-monotone curves */
  bool compare_curves(const X_monotone_curve_2& exp_answer,
                      const X_monotone_curve_2& real_answer)
  {
    typename Traits::Equal_2 equal = this->m_geom_traits.equal_2_object();
    if (equal(exp_answer, real_answer)) return true;

    std::string exp_answer_str = boost::lexical_cast<std::string>(exp_answer);
    std::string real_answer_str = boost::lexical_cast<std::string>(real_answer);
    this->print_answer(exp_answer_str, real_answer_str, "x-monotone curve");
    return false;
  }

  /*! Compare two numbers */
  template <typename Type_T>
  bool compare(const Type_T& exp_answer,
               const Type_T& real_answer,
               const char* str = "result")
  {
    if (exp_answer == real_answer) return true;
    std::string exp_answer_str = boost::lexical_cast<std::string>(exp_answer);
    std::string real_answer_str = boost::lexical_cast<std::string>(real_answer);
    this->print_answer(exp_answer_str, real_answer_str, str);
    return false;
  }

public:
  /*! Constructor */
  Traits_base_test(const Traits& traits);

  /*! Destructor */
  virtual ~Traits_base_test();

  /*! Parse the command line */
  virtual bool parse(int argc, char* argv[]);

  /*! Perform the test */
  virtual bool perform();

  /*! Clear the data structures */
  virtual void clear();
};

/*!
 * Constructor.
 * Accepts test data file name.
 */
template <typename Geom_traits_T>
Traits_base_test<Geom_traits_T>::Traits_base_test(const Traits& traits) :
  Base(traits),
  m_abort_on_error(false)       // run all tests
{
  m_violation_map[PRECONDITION] = std::string("precondition");
  m_violation_map[POSTCONDITION] = std::string("postcondition");
  m_violation_map[ASSERTION] = std::string("assertion");
  m_violation_map[WARNING] = std::string("warning");
}

/*! Destructor.
 */
template <typename Geom_traits_T>
Traits_base_test<Geom_traits_T>::~Traits_base_test() { clear(); }

template <typename Geom_traits_T>
bool Traits_base_test<Geom_traits_T>::parse(int argc, char* argv[])
{

  /* Waqar
  The arguments are
  argv 0 is ./test_traits (string)
  argv 1 is data/polycurves_conics/compare_y_at_x.pt
  argv 2 is data/polycurves_conics/compare_y_at_x.xcv
  argv 3 is data/polycurves_conics/compare_y_at_x.cv
  argv 4 is data/polycurves_conics/compare_y_at_x
  argv 5 is polycurve_conic_traits (string)
  */
  Base::parse(argc, argv);

  if (argc != 6) {
    this->print_info(std::string("Usage: ").append(argv[0]).
                     append(" points_file xcurves_file curves_file commands_file traits_type_name"));
    return false;
  }

  m_filename_commands.assign(argv[4]);
  m_traitstype.assign(argv[5]);

  return true;
}

/*! Clear the data structures */
template <typename Geom_traits_T>
void Traits_base_test<Geom_traits_T>::clear()
{
  Base::clear();
  m_filename_commands.clear();
}

/*!
 * Command dispatcher. Retrieves a line from the input file and performes
 * some action. See comments for suitable function in order to know specific
 * command arguments.
 */
template <typename Geom_traits_T>
bool Traits_base_test<Geom_traits_T>::perform()
{
  //std::cout << "*************" << m_filename_commands.c_str() << std::endl; // output is compare_y_at_x filepath
  std::ifstream is(m_filename_commands.c_str());

  if (!is.is_open())
  {
    this->print_error(std::string("cannot open file ").append(m_filename_commands));
    return false;
  }

  bool test_result = true;
  std::cout << "Performing test: traits type is " << m_traitstype //m_traitstype = polycurve_conic_traits, arg 6 in script
            << ", input files are "
            << this->m_filename_points << " "
            << this->m_filename_xcurves << " "
            << this->m_filename_curves << " "
            << m_filename_commands << std::endl << std::endl;;

  this->m_eol_printed = true;
  std::string line;
  char buff[1024];
  // bool abort = false;
  int counter = 0;

  while (this->skip_comments(is, line))
  {
    std::istringstream str_stream(line, std::istringstream::in);
    buff[0] = '\0';
    str_stream.getline(buff, 1024, ' ');
    std::string str_command(buff);    //buff is "compare_y_at_x" which is a string taken from the start of every line of the file compare_y_at_x
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
      CGAL_USE(location);
      continue;
#endif
    }

    try
    {
      bool result;
      bool ignore = exec(str_stream, str_command, result);
      if (ignore) continue;
      if ((m_violation_tested != NON) &&
          (m_violation_tested != m_violation_occurred))
      {
        //violation is expected but it did not occur
        result = false;
        // if (m_abort_on_error) abort = true;
      }

      this->print_result(result);
      test_result &= result;
    }

    catch (CGAL::Precondition_exception /* e */)
    {
      if (m_violation_tested != PRECONDITION)
      {
        test_result = false;
        // if (m_abort_on_error) abort = true;
      }
    }

    catch (CGAL::Postcondition_exception /* e */)
    {
      if (m_violation_tested != POSTCONDITION)
      {
        test_result = false;
        // if (m_abort_on_error) abort = true;
      }
    }

    catch (CGAL::Warning_exception /* e */)
    {
      if (m_violation_tested != WARNING)
      {
        test_result = false;
        // if (m_abort_on_error) abort = true;
      }
    }

    catch (CGAL::Assertion_exception /* e */)
    {
      if (m_violation_tested != ASSERTION)
      {
        test_result = false;
        // if (m_abort_on_error) abort = true;
      }
    }

  } //while (this->skip_comments(is, line)) loop

  is.close();
  return test_result;
}

/*!
 */
template <typename Geom_traits_T>
bool Traits_base_test<Geom_traits_T>::
translate_boolean(std::string& str_value)
{
  if (str_value == "TRUE" || str_value == "true") return true;
  return false;
}

/*!
 */
template <typename Geom_traits_T>
unsigned int
Traits_base_test<Geom_traits_T>::translate_enumerator(std::string& str_value)
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
template <typename Geom_traits_T>
std::pair<typename Traits_base_test<Geom_traits_T>::Enum_type, unsigned int>
Traits_base_test<Geom_traits_T>::translate_int_or_text(std::string& str_value)
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
template <typename Geom_traits_T>
bool Traits_base_test<Geom_traits_T>::
get_expected_boolean(std::istringstream& str_stream)
{
  char buff[1024];
  str_stream.getline( buff, 1024, '.');
  buff[str_stream.gcount()] = '\0';
  std::string str_expres = this->remove_blanks(buff);
  return translate_boolean(str_expres);
}

/*!
 */
template <typename Geom_traits_T>
unsigned int Traits_base_test<Geom_traits_T>::
get_expected_enum(std::istringstream& str_stream)
{
  char buff[1024];
  str_stream.getline(buff, 1024, '.');
  buff[str_stream.gcount()] = '\0';
  std::string str_expres = this->remove_blanks(buff);
  return translate_enumerator(str_expres);
}

/*!
 */
template <typename Geom_traits_T>
std::pair<typename Traits_base_test<Geom_traits_T>::Enum_type, unsigned int>
Traits_base_test<Geom_traits_T>::get_next_input(std::istringstream& str_stream)
{
  char buff[1024];
  do {
    str_stream.getline(buff, 1024, ' ');
  } while (str_stream.gcount() == 1);
  buff[str_stream.gcount()] = '\0';
  std::string str_expres = this->remove_blanks(buff);
  return translate_int_or_text(str_expres);
}

#endif
