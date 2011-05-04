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

template <class T_Traits>
class Traits_base_test {
protected:
  typedef T_Traits                                      Traits;
  typedef typename Traits::Point_2                      Point_2;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename Traits::Curve_2                      Curve_2;

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

  /*! An instance of the traits */
  Traits m_traits;

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
    typename Traits::Equal_2 equal = m_traits.equal_2_object();
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
    typename Traits::Equal_2 equal = m_traits.equal_2_object();
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

  template <class stream>
  bool read_point(stream& is, Point_2&);

  template <class stream>
  bool read_xcurve(stream& is, X_monotone_curve_2&);

  template <class stream>
  bool read_curve(stream& is, Curve_2&);

  //@{

  /*! Test Compare_x_2
   */
  bool compare_x_wrapper(std::istringstream& );

  /*! Test Compare_x_near_limit_2
   */
  bool compare_x_near_limit_wrapper(std::istringstream& );
  bool compare_x_near_limit_wrapper_imp(std::istringstream&,
                                        CGAL::Arr_use_dummy_tag);
  bool compare_x_near_limit_wrapper_imp(std::istringstream&,
                                        CGAL::Arr_use_traits_tag);
  
  /*! Test Compare_x_at_limit_2
   */
  bool compare_x_at_limit_wrapper(std::istringstream& );
  bool compare_x_at_limit_wrapper_imp(std::istringstream&,
                                      CGAL::Arr_use_dummy_tag);
  bool compare_x_at_limit_wrapper_imp(std::istringstream&,
                                      CGAL::Arr_use_traits_tag);

  /*! Test Compare_x_near_boundary_2
   */
  bool compare_x_near_boundary_wrapper(std::istringstream& );
  bool compare_x_near_boundary_wrapper_imp(std::istringstream&,
                                           CGAL::Arr_use_dummy_tag);
  bool compare_x_near_boundary_wrapper_imp(std::istringstream&,
                                           CGAL::Arr_use_traits_tag);
  
  /*! Test Compare_x_on_boundary_2
   */
  bool compare_x_on_boundary_wrapper(std::istringstream& );
  bool compare_x_on_boundary_wrapper_imp(std::istringstream&,
                                         CGAL::Arr_use_dummy_tag);
  bool compare_x_on_boundary_wrapper_imp(std::istringstream&,
                                         CGAL::Arr_use_traits_tag);
  
  /*! Test Compare_y_near_boundary_2
   */
  bool compare_y_near_boundary_wrapper(std::istringstream& );
  bool compare_y_near_boundary_wrapper_imp(std::istringstream&,
                                           CGAL::Arr_use_dummy_tag);
  bool compare_y_near_boundary_wrapper_imp(std::istringstream&,
                                           CGAL::Arr_use_traits_tag);

  /*! Test Parameter_space_in_x_2
   */
  bool parameter_space_in_x_wrapper(std::istringstream& );
  bool parameter_space_in_x_wrapper_imp(std::istringstream&,
                                        CGAL::Arr_use_dummy_tag);
  bool parameter_space_in_x_wrapper_imp(std::istringstream&,
                                        CGAL::Arr_use_traits_tag);
  
  /*! Test Parameter_space_in_y_2
   */
  bool parameter_space_in_y_wrapper(std::istringstream& );
  bool parameter_space_in_y_wrapper_imp(std::istringstream&,
                                        CGAL::Arr_use_dummy_tag);
  bool parameter_space_in_y_wrapper_imp(std::istringstream&,
                                        CGAL::Arr_use_traits_tag);

  /*! Compare_xy_2
   */
  bool compare_xy_wrapper(std::istringstream& line);

  /*! Tests Construct_min_vertex_2.
   *  Degenerate case: vertical curve.
   */
  bool min_vertex_wrapper(std::istringstream& line);

  /*! Tests Construct_max_vertex_2.
   * Degenerate case: vertical curve.
   */
  bool max_vertex_wrapper(std::istringstream& line);

  bool is_vertical_wrapper(std::istringstream& line);

  /*! Tests Compare_y_at_x_2.
   * Return the location of the given point with respect to the input curve.
   * Degenerate cases: The point is an endpoint of the curve.
   *                   The curve is vertical.
   */
  bool compare_y_at_x_wrapper(std::istringstream& line);

  /*! Tests Compare_y_at_x_left_2.
   * Compare the y value of two x-monotone curves immediately to the left
   * (resp. right) of their intersection point.
   * Degenerate cases: The curves coincide.
   *                   The curves coincide and vertical.
   *                   One of the curves is vertical.
   */
  bool compare_y_at_x_left_wrapper(std::istringstream& line);
  bool compare_y_at_x_left_wrapper_imp(std::istringstream& line,
                                       CGAL::Tag_false);
  bool compare_y_at_x_left_wrapper_imp(std::istringstream& line,
                                       CGAL::Tag_true);

  /*! Tests Compare_y_at_x_right_2.
   * Compare the y value of two x-monotone curves immediately to the right
   * (resp. right) of their intersection point.
   * Degenerate cases: The curves coincide.
   *                   The curves coincide and vertical.
   *                   One of the curves is vertical.
   */
  bool compare_y_at_x_right_wrapper(std::istringstream& line);

  /*! Tests Equal_2::operator()(Point_2, Point_2).
   * Check whether two points are the same.
   */
  bool equal_points_wrapper(std::istringstream& line);

  /*! Tests Equal_2::operator()(X_monotone_curve_2, X_monotone_curve_2).
   * Check whether two x-monotone curves are the same.
   */
  bool equal_curves_wrapper(std::istringstream& line);

  /*! Tests Make_x_monotone_2.
   * Cut the given curve into x-monotone subcurves and insert them into the
   * given output iterator.
   * Degenerate cases for polylines: The first segment is vertical. The last
   * segment is vertical. Both firt and last are vertical. An internal segment
   * is vertical.
   */
  bool make_x_monotone_wrapper(std::istringstream& line);

  /*! Tests Intersect_2.
   * Find the intersections of the two given curves and insert them into the 
   * given output iterator.
   * Degenerate cases for polylines: The most right (resp. left) endpoints of
   * the two curves coincide. Both endpoints coincide. The most right (resp.
   * left) endpoint of one curve and the first (resp. last) segment of the
   * other coincide.
   */
  bool intersect_wrapper(std::istringstream& line);

  /*! Tests Split_2.
   * Split a given x-monotone curve at a given point into two sub-curves.
   * Degenerate cases for polylines: the point and a polyline internal point
   * coincides.
   */
  bool split_wrapper(std::istringstream& line);

  /*! Tests Are_mergeable_2.
   * Check whether it is possible to merge two given x-monotone curves.
   */
  bool are_mergeable_wrapper(std::istringstream& line);
  bool are_mergeable_wrapper_imp(std::istringstream& line, CGAL::Tag_false);
  bool are_mergeable_wrapper_imp(std::istringstream& line, CGAL::Tag_true);

  /*! Tests Merge_2.
   * Merge two given x-monotone curves into a single curve.
   */
  bool merge_wrapper(std::istringstream& line);
  bool merge_wrapper_imp(std::istringstream& line, CGAL::Tag_false);
  bool merge_wrapper_imp(std::istringstream& line, CGAL::Tag_true);

  /*! Tests Approximate_2.
   * Return an approximation of a point coordinate.
   */
  bool approximate_wrapper(std::istringstream& line);

  /*! tests Construct_x_monotone_curve_2.
   * Return an x-monotone curve connecting the two given endpoints.
   */
  bool construct_x_monotone_curve_wrapper(std::istringstream& line);

  //@}

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

template <class T_Traits>
template <class stream>
bool Traits_base_test<T_Traits>::
read_point(stream& is, typename T_Traits::Point_2& p)
{
  Basic_number_type x, y;
  is >> x >> y;
  p = typename T_Traits::Point_2(x, y);
  return true;
}

template <class T_Traits>
template <class stream>
bool Traits_base_test<T_Traits>::
read_xcurve(stream& is, typename T_Traits::X_monotone_curve_2& xcv)
{
  Basic_number_type x1, y1, x2, y2;
  is >> x1 >> y1 >> x2 >> y2;
  Point_2 p1(x1, y1);
  Point_2 p2(x2, y2);
  CGAL_assertion(p1 != p2);
  xcv = typename T_Traits::X_monotone_curve_2(p1, p2);
  return true;
}

template <class T_Traits>
template <class stream>
bool
Traits_base_test<T_Traits>::read_curve(stream& is,
                                       typename T_Traits::Curve_2& cv)
{
  Basic_number_type x1, y1, x2, y2;
  is >> x1 >> y1 >> x2 >> y2;
  Point_2 p1(x1, y1);
  Point_2 p2(x2, y2);
  CGAL_assertion(p1 != p2);
  cv = typename T_Traits::Curve_2(p1, p2);
  return true;
}

#if TEST_TRAITS == CORE_CONIC_TRAITS || \
    TEST_TRAITS == RATIONAL_ARC_TRAITS

// conic traits and rational traits use same number 
// type CORE:Expr so this code can be shared

/*! Read a point */

template <>
template <class stream>
bool
Traits_base_test<Traits>::read_point(stream& is, Point_2& p)
{
  Rational rat_x,rat_y;
  is >> rat_x >> rat_y;
  Basic_number_type x(rat_x), y(rat_y);
  p = Point_2(x, y);
  return true;
}

#endif

#if TEST_TRAITS == SPHERICAL_ARC_TRAITS

/*! Read a point */

template <>
template <class stream>
bool
Traits_base_test<Traits>::read_point(stream& is, Point_2& p)
{
  Basic_number_type x, y, z;
  is >> x >> y >> z;
  p = Point_2(x, y, z);
  return true;
}

/*! Read a xcurve */
template <>
template <class stream>
bool
Traits_base_test<Traits>::read_xcurve(stream& is, X_monotone_curve_2& xcv)
{
  Point_2 p1,p2;
  read_point(is, p1);
  read_point(is, p2);
  CGAL_assertion(p1 != p2);
  xcv = X_monotone_curve_2(p1, p2);
  return true;
}

/*! Read a curve */
template <>
template <class stream>
bool
Traits_base_test<Traits>::read_curve(stream& is, Curve_2& cv)
{
  Point_2 p1, p2;
  read_point(is, p1);
  read_point(is, p2);
  CGAL_assertion(p1 != p2);
  cv = Curve_2(p1, p2);
  return true;
}

#endif

#if TEST_TRAITS == RATIONAL_ARC_TRAITS

/*! Read a xcurve */
template <>
template <class stream>
bool
Traits_base_test<Traits>::read_xcurve(stream& is, X_monotone_curve_2& xcv)
{
  Curve_2 tmp_cv;
  if (!read_curve(is, tmp_cv))
    return false;
  xcv = X_monotone_curve_2(tmp_cv);
  return true;
}

template <class stream>
bool read_coefficients(stream& is, Rat_vector& coeffs)
{
  unsigned int num_coeffs;
  Rational rat;
  is >> num_coeffs;
  coeffs.clear();
  for (unsigned int j = 0; j < num_coeffs; j++) {
    is >> rat;
    coeffs.push_back(rat);
  }
  return true;
}

/*! Read a curve */
template <>
template <class stream>
bool Traits_base_test<Traits>::read_curve(stream& is, Curve_2& cv)
{
  // Get the arc type:
  Rat_vector p_coeffs, q_coeffs;
  Algebraic src, trg;
  int dir = 0;
  char type;
  is >> type;
  if (type == 'a' || type == 'A') 
  {
    //Default constructor
    cv = Curve_2();
    return true;
  }
  else if (type == 'b' || type == 'B') 
  {
    //Constructor of a whole polynomial curve
    if (read_coefficients(is,p_coeffs))
      cv = Curve_2(p_coeffs);
    else
      return false;
    return true;
  }
  else if (type == 'c' || type == 'C') 
  {
    //Constructor of a polynomial ray
    if (!read_coefficients(is,p_coeffs))
      return false;
    is >> src >> dir;
    cv = Curve_2(p_coeffs, src, (dir == 0 ? false : true));
    return true;
  }
  else if (type == 'd' || type == 'D') 
  {
    //Constructor of a polynomial arc
    if (!read_coefficients(is,p_coeffs))
      return false;
    is >> src >> trg;
    cv = Curve_2(p_coeffs, src, trg);
    return true;
  }
  else if (type == 'e' || type == 'E') 
  {
    //Constructor of a whole rational function
    if (!read_coefficients(is,p_coeffs))
      return false;
    if (!read_coefficients(is,q_coeffs))
      return false;
    cv = Curve_2(p_coeffs, q_coeffs);
    return true;
  }
  else if (type == 'f' || type == 'F') 
  {
    //Constructor of a ray of a rational function
    if (!read_coefficients(is,p_coeffs))
      return false;
    if (!read_coefficients(is,q_coeffs))
      return false;
    is >> src >> dir;
    cv = Curve_2(p_coeffs, q_coeffs, src, (dir == 0 ? false : true));
    return true;
  }
  else if (type == 'g' || type == 'G') 
  {
    //Constructor of a bounded rational arc
    if (!read_coefficients(is, p_coeffs))
      return false;
    if (!read_coefficients(is, q_coeffs))
      return false;
    is >> src >> trg;
    cv = Curve_2(p_coeffs, q_coeffs, src, trg);
    return true;
  }
  // If we reached here, we have an unknown rational arc type:
  std::cerr << "Illegal rational arc type specification: " << type << "."
            << std::endl;
  return (false);
}

#endif

#if TEST_TRAITS == POLYLINE_TRAITS || TEST_TRAITS == NON_CACHING_POLYLINE_TRAITS

template <>
template <class stream>
bool
Traits_base_test<Traits>::
read_xcurve(stream& is, Traits::X_monotone_curve_2& xcv)
{
  unsigned int num_points;
  is >> num_points;
  std::vector<Point_2> points;
  points.clear();
  for (unsigned int j = 0; j < num_points; j++) {
    Basic_number_type x, y;
    is >> x >> y;
    Point_2 p(x, y);
    points.push_back(p);
  }
  xcv =
    CGAL::
    Arr_polyline_traits_2<Segment_traits>::X_monotone_curve_2(points.begin(),
                                                              points.end());
  return true;
}

template <>
template <class stream>
bool
Traits_base_test<Traits>::read_curve(stream& is, Traits::Curve_2& cv)
{
  unsigned int num_points;
  is >> num_points;
  std::vector<Point_2> points;
  points.clear();
  for (unsigned int j = 0; j < num_points; j++) {
    Basic_number_type x, y;
    is >> x >> y;
    Point_2 p(x, y);
    points.push_back(p);
  }
  cv = CGAL::Arr_polyline_traits_2<Segment_traits>::Curve_2(points.begin(),
                                                            points.end());
  return true;
}

#elif TEST_TRAITS == LINEAR_TRAITS

template <>
template <class stream>
bool
Traits_base_test<Traits>::read_xcurve(stream& is, X_monotone_curve_2& xcv)
{
  is >> xcv;
  return true;
}

template <>
template <class stream>
bool
Traits_base_test<Traits>::read_curve(stream& is, Curve_2& cv)
{
  is >> cv;
  return true;
}

#elif TEST_TRAITS == CORE_CONIC_TRAITS

/*! */
template <class stream>
bool read_orientation(stream& is, CGAL::Orientation& orient)
{
  int i_orient;
  is >> i_orient;
  orient = (i_orient > 0) ? CGAL::COUNTERCLOCKWISE :
    (i_orient < 0) ? CGAL::CLOCKWISE : CGAL::COLLINEAR;
  return true;
}

/*! */
template <class stream>
bool read_app_point(stream& is, Point_2& p)
{
  double x, y;
  is >> x >> y;
  p = Point_2(Algebraic(x), Algebraic(y));
  return true;
}

/*! */
template <class stream>
bool read_orientation_and_end_points(stream& is, CGAL::Orientation& orient,
                                     Point_2& source, Point_2& target)
{
  // Read the orientation.
  if (!read_orientation(is, orient)) return false;

  // Read the end points of the arc and create it.
  if (!read_app_point(is, source)) return false;
  if (!read_app_point(is, target)) return false;
  return true;
}

/*! Read a circle or an ellipse */
template <class stream>
bool read_ellipse(stream& is, bool& is_circle, Rat_circle& circle,
                  Rational& r, Rational& s,
                  Rational& t, Rational& u, Rational& v, Rational& w)
{
  // Read the ellipse (using the format "a b x0 y0"):
  //              2               2
  //   ( x - x0 )      ( y - y0 )
  //    --------   +   --------   = 1
  //       a              b
  //
  Rational a, b, x0, y0;
  is >> a >> b >> x0 >> y0;

  if (a == b) {
    is_circle = true;
    circle = Rat_circle (Rat_point (x0, y0), a);
  }
  else {
    r = 1/a;
    s = 1/b;
    t = 0;
    u = -2*x0*b;
    v = -2*y0*a;
    w = x0*x0*b + y0*y0*a - a*b;
  }
  return true;
}

/*! */
template <class stream, class Curve>
bool read_partial_ellipse(stream& is, Curve& cv)
{
  bool is_circle;               // Is this a circle.
  Rat_circle circle;
  Rational r, s, t, u, v, w;
  if (!read_ellipse(is, is_circle, circle, r, s, t, u, v, w)) return false;
  CGAL::Orientation orient;
  Point_2 source, target;
  if (!read_orientation_and_end_points(is, orient, source, target))  //!!!
    return false;

  // Create the conic (or circular) arc.
  cv = (is_circle) ? Curve (circle, orient, source, target) :
    Curve (r, s, t, u, v, w, orient, source, target);
  return true;
}

/*! */
template <class stream, class Curve>
bool read_full_ellipse(stream& is, Curve& cv)
{
  bool is_circle;               // Is this a circle.
  Rat_circle circle;
  Rational r, s, t, u, v, w;
  if (!read_ellipse(is, is_circle, circle, r, s, t, u, v, w))
    return false;

  // Create a full ellipse (or circle).
  cv = (is_circle) ? Curve (circle) : Curve (r, s, t, u, v, w);
  return true;
}

/*! Read a hyperbola */
template <class stream, class Curve>
bool read_hyperbola(stream& is, Curve& cv)
{
  // Read the hyperbola (using the format "a b x0 y0"):
  //              2              2
  //    ( x - x0 )     ( y - y0 )
  //     --------   -   --------   = 1
  //       a               b
  //
  Rational a, b, x0, y0;
  is >> a >> b >> x0 >> y0;

  Rational r = b;
  Rational s= -a;
  Rational t = 0;
  Rational u = -2*x0*b;
  Rational v = 2*y0*a;
  Rational w = x0*x0*b - y0*y0*a - a*b;

  CGAL::Orientation orient;
  Point_2 source, target;
  if (!read_orientation_and_end_points(is, orient, source, target)) //!!!
    return false;

  // Create the conic (or circular) arc.
  cv = Curve (r, s, t, u, v, w, orient, source, target);
  return true;
}

/*! Read a hyperbola */
template <class stream, class Curve>
bool read_parabola(stream& is, Curve& cv)
{
  // Read the parabola (using the format "c x0 y0"):
  //
  //          2
  //  (x - x0) = 4c*(y - y0)
  //
  Rational c, x0, y0;
  is >> c >> x0 >> y0;

  Rational r = 1;
  Rational s = 0;
  Rational t = 0;
  Rational u = -2*x0;
  Rational v = -4*c;
  Rational w = x0*x0 + 4*c*y0;

  CGAL::Orientation orient;
  Point_2 source, target;
  if (!read_orientation_and_end_points(is, orient, source, target))
    return false;

  // Create the conic (or circular) arc.
  cv = Curve (r, s, t, u, v, w, orient, source, target);
  return true;
}

/*! */
template <class stream, class Curve>
bool read_segment(stream& is, Curve& cv)
{

  // Read a segment, given by its endpoints (x1,y1) and (x2,y2);
  Rational x1, y1, x2, y2;
  is >> x1 >> y1 >> x2 >> y2;
  
  Rat_point source (x1, y1);
  Rat_point target (x2, y2);
  Rat_segment segment (source, target);

  // Create the segment.
  cv = Curve (segment);
  return true;
}

/*! */
template <class stream, class Curve>
bool read_general_arc(stream& is, Curve& cv)
{
  // Read a general conic, given by its coefficients <r,s,t,u,v,w>.
  Rational r, s, t, u, v, w;                // The conic coefficients.
  is >> r >> s >> t >> u >> v >> w;
    // Read the orientation.
  int i_orient = 0;
  is >> i_orient;
  CGAL::Orientation orient = (i_orient > 0) ? CGAL::COUNTERCLOCKWISE :
    (i_orient < 0) ? CGAL::CLOCKWISE : CGAL::COLLINEAR;

  // Read the approximated source, along with a general conic 
  // <r_1,s_1,t_1,u_1,v_1,w_1> whose intersection with <r,s,t,u,v,w>
  // defines the source.
  Point_2 app_source;
  if (!read_app_point(is, app_source)) return false;
  Rational r1, s1, t1, u1, v1, w1;
  is >> r1 >> s1 >> t1 >> u1 >> v1 >> w1;

  // Read the approximated target, along with a general conic 
  // <r_2,s_2,t_2,u_2,v_2,w_2> whose intersection with <r,s,t,u,v,w>
  // defines the target.
  Point_2 app_target;
  if (!read_app_point(is, app_target)) return false;

  Rational r2, s2, t2, u2, v2, w2;
  is >> r2 >> s2 >> t2 >> u2 >> v2 >> w2;

  // Create the conic arc.
  cv = Curve (r, s, t, u, v, w, orient,
              app_source, r1, s1, t1, u1, v1, w1,
              app_target, r2, s2, t2, u2, v2, w2);
  return true;
}

/*! */
template <class stream, class Curve>
bool read_general_curve(stream& is, Curve& cv)
{
  Rational r, s, t, u, v, w;                // The conic coefficients.
  // Read a general conic, given by its coefficients <r,s,t,u,v,w>.
  is >> r >> s >> t >> u >> v >> w;
  CGAL::Orientation orient;
  Point_2 source, target;
  if (!read_orientation_and_end_points(is, orient, source, target))
    return false;

  // Create the conic (or circular) arc.
  cv = Curve (r, s, t, u, v, w, orient, source, target);
  return true;
}

/*! Read an x-monotone conic curve */
template <>
template <class stream>
bool
Traits_base_test<Traits>::
read_xcurve(stream& is, X_monotone_curve_2& xcv)
{
  Curve_2 tmp_cv;
  if (!read_curve(is,tmp_cv))
    return false;
  xcv=X_monotone_curve_2(tmp_cv);
  return true;
}

/*! Read a general conic curve */
template <>
template <class stream>
bool
Traits_base_test<Traits>::read_curve(stream& is, Curve_2& cv)
{
  // Get the arc type:
  char type;
  is >> type;
  if (type == 'f' || type == 'F') 
  {
    return read_full_ellipse(is, cv);
  }
  else if (type == 's' || type == 'S') 
  {
    return read_segment(is, cv);
  }
  else if (type == 'i' || type == 'I') 
  {
    return read_general_arc(is, cv);
  }
  else if (type == 'c' || type == 'C') 
  {
    // Read a general conic, given by its coefficients <r,s,t,u,v,w>.
    Rational r, s, t, u, v, w;
    is >> r >> s >> t >> u >> v >> w;
    // Create a full conic (should work only for ellipses).
    cv = Curve_2(r, s, t, u, v, w);
    return (true);
  }
  else if (type == 'e' || type == 'E') 
  {
    return read_partial_ellipse(is, cv);
  }
  else if (type == 'h' || type == 'H') 
  {
    return read_hyperbola(is, cv);
  }
  else if (type == 'p' || type == 'P') 
  {
    return read_parabola(is, cv);
  }
  else if (type == 'a' || type == 'A') 
  {
    return read_general_curve(is, cv);
  }

  // If we reached here, we have an unknown conic type:
  std::cerr << "Illegal conic type specification: " << type << "."
	    << std::endl;
  return (false);
}

#elif TEST_TRAITS == CIRCLE_SEGMENT_TRAITS

template <class stream>
bool read_ort_point(stream& is, Point_2& p)
{
  bool is_rat;
  typename Point_2::CoordNT ort_x, ort_y;
  Number_type alpha,beta,gamma;
  is >> is_rat;
  if (is_rat)
  {
    is >> alpha;
    ort_x=Point_2::CoordNT(alpha);
  }
  else
  {
    is >> alpha >> beta >> gamma;
    ort_x=Point_2::CoordNT(alpha,beta,gamma);
  }
  is >> is_rat;
  if (is_rat)
  {
    is >> alpha;
    ort_y=Point_2::CoordNT(alpha);
  }
  else
  {
    is >> alpha >> beta >> gamma;
    ort_y=Point_2::CoordNT(alpha,beta,gamma);
  }
  p = Point_2(ort_x, ort_y);
  return true;
}

/*! Read an x-monotone circle segment curve */
template <>
template <class stream>
bool
Traits_base_test<Traits>::read_xcurve(stream& is,X_monotone_curve_2& xcv)
{
  bool ans=true;
  char type;
  is >> type;
  if (type == 'z' || type == 'Z')
  {
    Line_2 l;
    Point_2 ps,pt;
    is >> l;
    ans &= read_ort_point(is, ps);
    ans &= read_ort_point(is, pt);
    xcv=X_monotone_curve_2(l, ps, pt);
    return ans;
  }
  else if (type == 'y' || type == 'Y')
  {
    Rat_point_2 ps,pt;
    is >> ps >> pt;
    xcv=X_monotone_curve_2(ps,pt);
    return true;
  }
  else if (type == 'x' || type == 'X')
  {
    Circle_2 c;
    Point_2 ps,pt;
    is >> c;
    ans &= read_ort_point(is, ps);
    ans &= read_ort_point(is, pt);
    xcv=X_monotone_curve_2(c, ps, pt, c.orientation());
    return ans;
  }
  // If we reached here, we have an unknown conic type:
  std::cerr << "Illegal circle segment type specification: " << type
            << std::endl;
  return false;
}

/*! Read a general circle segment curve */
template <>
template <class stream>
bool
Traits_base_test<Traits>::read_curve(stream& is,Curve_2& cv)
{
  bool ans=true;
  char type;
  is >> type;
  if (type == 'a' || type == 'A')
  {
    Rat_point_2 ps,pt;
    is >> ps >> pt;
    Segment_2 s(ps,pt);
    cv=Curve_2(s);
    return true;
  }
  else if (type == 'b' || type == 'B')
  {
    Rat_point_2 ps,pt;
    is >> ps >> pt;
    cv=Curve_2(ps,pt);
    return true;
  }
  else if (type == 'c' || type == 'C')
  {
    Line_2 l;
    Point_2 ps,pt;
    is >> l;
    ans &= read_ort_point(is, ps);
    ans &= read_ort_point(is, pt);
    cv=Curve_2(l,ps,pt);
    return ans;
  }
  else if (type == 'd' || type == 'D')
  {
    Circle_2 c;
    is >> c;
    cv=Curve_2(c);
    return true;
  }
  else if (type == 'e' || type == 'E')
  {
    Rat_point_2 p;
    Rat_nt r;
    int orient;
    is >> p >> r >> orient;
    cv=Curve_2(p,r,static_cast<CGAL::Orientation>(orient));
    return true;
  }
  else if (type == 'f' || type == 'F')
  {
    Circle_2 c;
    Point_2 ps,pt;
    is >> c;
    ans &= read_ort_point(is, ps);
    ans &= read_ort_point(is, pt);
    cv=Curve_2(c,ps,pt);
    return ans;
  }
  else if (type == 'g' || type == 'G')
  {
    Rat_point_2 p;
    Rat_nt r;
    int orient;
    Point_2 ps,pt;
    is >> p >> r >> orient;
    ans &= read_ort_point(is, ps);
    ans &= read_ort_point(is, pt);
    cv=Curve_2(p,r,static_cast<CGAL::Orientation>(orient),ps,pt);
    return ans;
  }
  else if (type == 'h' || type == 'H')
  {
    Rat_point_2 ps,pm,pt;
    is >> ps >> pm >> pt;
    cv=Curve_2(ps,pm,pt);
    return true;
  }
  // If we reached here, we have an unknown conic type:
  std::cerr << "Illegal circle segment type specification: " << type
            << std::endl;
  return false;
}

#elif TEST_TRAITS == BEZIER_TRAITS

template <>
template <class stream>
bool
Traits_base_test<Traits>::read_point(stream& is, Point_2& p)
{
  Rational rat_x,rat_y;
  is >> rat_x >> rat_y;
  p = Point_2(rat_x, rat_y);
  return true;
}

/*! Read an x-monotone bezier curve */

template <>
template <class stream>
bool
Traits_base_test<Traits>::read_xcurve(stream& is, X_monotone_curve_2& xcv)
{
  std::list<CGAL::Object>                  x_objs;
  std::list<CGAL::Object>::const_iterator  xoit;
  Curve_2 tmp_cv;
  is >> tmp_cv;
  Rational B_psx = Rational(tmp_cv.control_point(0).x());
  Rational B_psy = Rational(tmp_cv.control_point(0).y());
  Rational B_ptx =
    Rational(tmp_cv.control_point(tmp_cv.number_of_control_points()-1).x());
  Rational B_pty =
    Rational(tmp_cv.control_point(tmp_cv.number_of_control_points()-1).y());
  Point_2 B_ps(B_psx, B_psy);
  Point_2 B_pt(B_ptx, B_pty);
  Traits::Make_x_monotone_2 make_x_monotone =
    this->m_traits.make_x_monotone_2_object();
  make_x_monotone (tmp_cv, std::front_inserter (x_objs));
  xoit = x_objs.begin();
  if (CGAL::assign (xcv, *xoit))
    return true;
  return false;
}

/*! Read a general bezier curve */
template <>
template <class stream>
bool
Traits_base_test<Traits>::read_curve(stream& is, Curve_2& cv)
{
  is >> cv;
  return true;
}

#elif TEST_TRAITS == LINE_ARC_TRAITS || \
  TEST_TRAITS == CIRCULAR_ARC_TRAITS || \
  TEST_TRAITS == CIRCULAR_LINE_ARC_TRAITS

/*! Read an arc point */
template <typename T_Traits, typename stream>
bool read_arc_point(stream& is, typename T_Traits::Point_2& p)
{
  Basic_number_type x, y;
  is >> x >> y;
  Circular_kernel::Point_2 lp(x, y);
  p = typename T_Traits::Point_2(lp);
  return true;
}

bool is_deg_1(char c)
{
  return (c=='z' || c=='Z') || (c=='y' || c=='Y') || (c=='x' || c=='X') ||
         (c=='w' || c=='W') || (c=='v' || c=='V') || (c=='l' || c=='L');
}

bool is_deg_2(char c)
{
  return (c=='b' || c=='B') || (c=='c' || c=='C') ||
         (c=='d' || c=='D') || (c=='e' || c=='E');
}

#if TEST_TRAITS == LINE_ARC_TRAITS || \
  TEST_TRAITS == CIRCULAR_LINE_ARC_TRAITS

template <class stream>
Circular_kernel::Line_arc_2 read_line(char type, stream& is)
{
  if (type == 'z' || type == 'Z')
  {
    Circular_kernel::Line_2 l_temp;
    Circular_kernel::Circle_2 c_temp1,c_temp2;
    bool b1,b2;
    is >> l_temp >> c_temp1 >> b1 >> c_temp2 >> b2;
    return Circular_kernel::Line_arc_2(l_temp,c_temp1,b1,c_temp2,b2);
  }
  else if (type == 'y' || type == 'Y')
  {
    Circular_kernel::Line_2 l_temp,l_temp1,l_temp2;
    is >> l_temp >> l_temp1 >> l_temp2;
    return Circular_kernel::Line_arc_2(l_temp,l_temp1,l_temp2);
  }
  else if (type == 'x' || type == 'X')
  {
    Circular_kernel::Line_2 l_temp;
    Circular_kernel::Circular_arc_point_2 p0,p1;
    is >> l_temp >> p0 >> p1;
    //std::cout << "got here l_temp p0 p1 " << l_temp << " " << p0 << " " << p1 << std::endl;
    return Circular_kernel::Line_arc_2(l_temp, p0, p1);
  }
  else if (type == 'w' || type == 'W' || type == 'l' || type == 'L')
  {
    Circular_kernel::Point_2 p0,p1;
    is >> p0 >> p1;
    return Circular_kernel::Line_arc_2(p0,p1);
  }
  else if (type == 'v' || type == 'V')
  {
    Circular_kernel::Segment_2 seg;
    is >> seg;
    return Circular_kernel::Line_arc_2(seg);
  }
  std::cout << "should never happen Line_arc_2 " << type <<std::endl;
  return Circular_kernel::Line_arc_2(); //should never happen
}
#endif

#if TEST_TRAITS == CIRCULAR_ARC_TRAITS || \
  TEST_TRAITS == CIRCULAR_LINE_ARC_TRAITS
template <class stream>
Circular_kernel::Circular_arc_2 read_arc(char type,stream& is)
{
  if (type == 'b' || type == 'B')
  {
    Circular_kernel::Circle_2 circle, circle1, circle2;
    bool b1, b2;
    is >> circle >> circle1 >> b1 >> circle2 >> b2;

    return Circular_kernel::Circular_arc_2(circle, circle1, b1, circle2, b2);
  }
  else if (type == 'c' || type == 'C')
  {
    Circular_kernel::Circle_2 circle;
    Circular_kernel::Circular_arc_point_2 p0, p1;
    is >> circle >> p0 >> p1;
    return Circular_kernel::Circular_arc_2(circle, p0, p1);
  }
  else if (type == 'd' || type == 'D')
  {
    Circular_kernel::Circle_2 circle;
    Circular_kernel::Line_2 line1, line2;
    bool b1,b2;
    is >> circle >> line1 >> b1 >> line2 >> b2;
    return Circular_kernel::Circular_arc_2(circle, line1, b1, line2, b2);
  }
  else
    CGAL_error_msg("Unrecognized constructor. Should never happen"      \
                   "Circular_arc_2");
  //  else if (type == 'e' || type == 'E')
  //  {
  //    Circular_kernel::Circular_arc_2 arc;
  //    Circular_kernel::Circle_2 circle;
  //    bool b1, b2;
  //    is >> arc >> b1 >> circle >> b2;
  //    return Circular_kernel::Circular_arc_2(arc, b1, circle, b2);
  //  }
  return Circular_kernel::Circular_arc_2(); //should never happen
}
#endif

#if TEST_TRAITS == LINE_ARC_TRAITS

/*! Read a line arc point */
template <>
template <class stream>
bool
Traits_base_test<Traits>::read_point(stream& is, Point_2& p)
{
  return read_arc_point<Traits, stream>(is, p);
}

/*! Read an x-monotone line arc curve */
template <>
template <class stream>
bool
Traits_base_test<Traits>::read_xcurve(stream& is, X_monotone_curve_2& xcv)
{
  // Get the arc type:
  char type;
  is >> type;
  if (is_deg_1(type))
  {
    xcv=read_line(type,is);
    return true;
  }
  return false;

}

/*! Read a general line arc curve */
template <>
template <class stream>
bool
Traits_base_test<Traits>::read_curve(stream& is, Curve_2& cv)
{
  // Get the arc type:
  char type;
  is >> type;
  if (is_deg_1(type))
  {
    cv=read_line(type,is);
    return true;
  }
  return false;
}

#elif TEST_TRAITS == CIRCULAR_ARC_TRAITS

/*! Read a circular arc point */
template <>
template <class stream>
bool
Traits_base_test<Traits>::read_point(stream& is, Point_2& p)
{
  return read_arc_point<Traits, stream>(is, p);
}

/*! Read an x-monotone circular arc curve */
template <>
template <class stream>
bool
Traits_base_test<Traits>::read_xcurve(stream& is,X_monotone_curve_2& xcv)
{
  // Get the arc type:
  char type;
  is >> type;
  if (is_deg_2(type))
  {
    xcv=read_arc(type, is);
    return true;
  }
  return false;
}

/*! Read a general circular curve */
template <>
template <class stream>
bool
Traits_base_test<Traits>::read_curve(stream& is, Curve_2& cv)
{
  // Get the arc type:
  char type;
  is >> type;
  if (type == 'a' || type == 'A') 
  {
    Circular_kernel::Circle_2 circle;
    is >> circle;
    cv=Circular_kernel::Circular_arc_2(circle);
    return true;
  }
  else if (is_deg_2(type))
  {
    cv=read_arc(type, is);
    return true;
  }
  return false;
}

#elif TEST_TRAITS == CIRCULAR_LINE_ARC_TRAITS

/*! Read a circular-line arc point */
template <>
template <class stream>
bool
Traits_base_test<Traits>::read_point(stream& is, Point_2& p)
{
  return read_arc_point<Traits, stream>(is, p);
}

/*! Read an x-monotone circular-line arc curve */
template <>
template <class stream>
bool
Traits_base_test<Traits>::read_xcurve(stream& is, X_monotone_curve_2& xcv)
{
  // Get the arc type:
  char type;
  is >> type;
  if (is_deg_1(type)) {
    xcv = read_line(type,is);
    return true;
  }
  else if (is_deg_2(type)) {
    xcv = X_monotone_curve_2(read_arc(type, is));
    return true;
  }
  return false;
}

/*! Read a general circular-line curve */
template <>
template <class stream>
bool Traits_base_test<Traits>::read_curve(stream& is, Curve_2& cv)
{
  // Get the arc type:
  char type;
  is >> type;
  if (type == 'a' || type == 'A') {
    Circular_kernel::Circle_2 circle;
    is >> circle;
    cv=Curve_2(circle);
    return true;
  }
  else if (is_deg_1(type)) {
    cv = Curve_2(read_line(type, is));
    return true;
  }
  else if (is_deg_2(type)) {
    cv = read_arc(type, is);
    return true;
  }
  return false;

}

#endif

#endif

#if TEST_TRAITS == ALGEBRAIC_TRAITS

#include <CGAL/IO/io.h>

template <>
template <class stream>
bool
Traits_base_test<Traits>::read_point(stream& is, Point_2& p) {
  Traits traits;
  Traits::Construct_point_2 construct_point_2
    = traits.construct_point_2_object();
  char type;
  is >> type;
  switch(type) {
  case 'i': {
    int x=0,y=0;
    is >> x >> y;
    p=construct_point_2(x, y);
    break;
  }
  case 'b': {
    Traits::Bound x,y;
    is >> x >> y;
    p=construct_point_2(x, y);
    break;
  }
  case 'c': {
    Traits::Coefficient x, y;
    is >> x >> y;
    p=construct_point_2(x, y);
    break;
  }
  case 'a': {
    Traits::Algebraic_real_1 x, y;
    is >> x >> y;
    p=construct_point_2(x, y);
    break;
  }
  case 's': {
    Traits::Algebraic_real_1 x;
    is >> x;
    Traits::X_monotone_curve_2 xcv;
    CGAL::swallow(is,'(');
    CGAL_assertion_code(bool check=)
      read_xcurve(is, xcv);
    CGAL_assertion(check);
    
    CGAL::swallow(is,')');
    p=construct_point_2(x, xcv);
    break;
  }
  case 'g': {
    Traits::Algebraic_real_1 x;
    is >> x;
    Traits::Curve_2 c;
    CGAL_assertion_code(bool check = )
      read_curve(is,c);
    CGAL_assertion(check);
    int arcno=0;
    is >> arcno;
    p=construct_point_2(x, c, arcno);
    break;
  }
  default: {
    std::cout << "Expected i, b, c, a, s, or g, but got \"" << type << "\""
	      << std::endl;
    return false;
  }
  }
  return true;
}

template <>
template <class stream>
bool Traits_base_test<Traits>::read_xcurve(stream& is, 
                                           Traits::X_monotone_curve_2& xcv)
{
  Traits traits;
  Traits::Construct_x_monotone_segment_2 construct_segment_2
    = traits.construct_x_monotone_segment_2_object();
  char type;
  is >> type;
  switch(type) {
   case '1': {
     Curve_2 cv;
     Point_2 end_left,end_right;
     CGAL_assertion_code(bool check=)
     read_curve(is,cv);
     CGAL_assertion(check);
     CGAL::swallow(is,'(');
     CGAL_assertion_code(check=)
     read_point(is,end_left);
     CGAL_assertion(check);
     CGAL::swallow(is,')');
     CGAL::swallow(is,'(');
     CGAL_assertion_code(check=)
     read_point(is,end_right);
     CGAL_assertion(check);
     CGAL::swallow(is,')');
     std::vector<Traits::X_monotone_curve_2> xcvs;
     construct_segment_2(cv, end_left, end_right, std::back_inserter(xcvs));
     CGAL_assertion(xcvs.size() == 1);
     xcv = xcvs[0];
     break;
    }
   case '2': {
     Curve_2 cv;
     Point_2 p;
     CGAL_assertion_code(bool check=)
     read_curve(is,cv);
     CGAL_assertion(check);
     CGAL::swallow(is,'(');
     CGAL_assertion_code(check=)
     read_point(is,p);
     CGAL_assertion(check);
     CGAL::swallow(is,')');
     std::string site_of_p_string;
     Traits::Site_of_point site_of_p;
     is >> site_of_p_string;
     if (site_of_p_string=="MIN_ENDPOINT") {
       site_of_p=Traits::MIN_ENDPOINT;
     } else if (site_of_p_string=="MAX_ENDPOINT") {
       site_of_p=Traits::MAX_ENDPOINT;
     } else {
       CGAL_assertion(site_of_p_string=="POINT_IN_INTERIOR");
       site_of_p=Traits::POINT_IN_INTERIOR;
     }
     std::vector<Traits::X_monotone_curve_2> xcvs;
     construct_segment_2(cv, p, site_of_p, std::back_inserter(xcvs));
     CGAL_assertion(xcvs.size() == 1);
     xcv = xcvs[0];
     break;
    }
   default: {
     std::cout << "Expected 1 or 2, but got \"" << type << "\"" << std::endl;
     return false;
   }
  }
  return true;
}

template <>
template <class stream>
bool Traits_base_test<Traits>::read_curve(stream& is, Curve_2& cv) {
  Traits traits;
  Traits::Polynomial_2 p;
  Traits::Construct_curve_2 construct_curve_2 
    = traits.construct_curve_2_object();
  is >> p;
  cv = construct_curve_2(p);
  return true;
}

#endif

#endif
