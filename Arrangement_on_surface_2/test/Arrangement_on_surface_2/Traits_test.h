#ifndef CGAL_TRAITS_TEST_H
#define CGAL_TRAITS_TEST_H

#include <CGAL/Object.h>
#include <CGAL/tags.h>
#include <CGAL/exceptions.h>

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>

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

enum Enum_type {NUMBER, SIGN, CURVE_END, BOUNDARY};

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
    Test_exception(Exception_type e_tmp) :
      CGAL::Failure_exception(global_lib, global_expr, global_file,
                              global_line, global_msg, global_type)
    {
      e=e_tmp;
    }
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
  std::string tmp=std::string(file);
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
  global_type=std::string(( !type ? "" : type ));
  global_expr=std::string(( !expr ? "" : expr ));
  global_line=line;
  global_msg=std::string(( !msg ? "" : msg ));
  if (!global_type.compare("precondition"))
  {
    violation_occurred=PRECONDITION;
    throw_exceptions(violation_tested==PRECONDITION);
  }
  else if (!global_type.compare("postcondition"))
  {
    violation_occurred=POSTCONDITION;
    throw_exceptions(violation_tested==POSTCONDITION);
  }
  else if (!global_type.compare("assertion"))
  {
    violation_occurred=ASSERTION;
    throw_exceptions(violation_tested==ASSERTION);
  }
  else if (!global_type.compare("warning"))
  {
    violation_occurred=WARNING;
    throw_exceptions(violation_tested==WARNING);
  }
}

CGAL::Failure_function prev_error_handler;
CGAL::Failure_function prev_warning_handler;

template <class T_Traits>
class Traits_test {
private:
  typedef T_Traits                                      Traits;
  typedef typename Traits::Point_2                      Point_2;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename Traits::Curve_2                      Curve_2;

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

  /*! A map between (strings) commands and (member functions) operations */
  typedef bool (Traits_test::* Wrapper)(std::istringstream &);
  typedef std::map<std::string, Wrapper>        Wrapper_map;
  typedef typename Wrapper_map::iterator        Wrapper_iter;
  Wrapper_map m_wrappers;

  /*! Collect the data depending on obj_t */
  bool read_input(std::ifstream & is , Read_object_type obj_t);

  /*! Perform the test */
  bool perform(std::ifstream & is);

  /*! Skip comments */
  void skip_comments(std::ifstream & is, char * one_line);

  std::string remove_blanks(char * str);
  bool get_expected_boolean(std::istringstream & str_stream);
  unsigned int get_expected_enum(std::istringstream & str_stream);
  bool translate_boolean(std::string & str_value);
  unsigned int translate_enumerator(std::string & str_value);
  std::pair<Enum_type,unsigned int> translate_int_or_text(std::string & str_value);//Enum_type
  std::pair<Enum_type,unsigned int> get_next_input(std::istringstream & str_stream);//Enum_type

  bool compare_points(const Point_2 & exp_answer, const Point_2 & real_answer,
               const char * str = "result")
  {
    typename Traits::Equal_2 equal = m_traits.equal_2_object();
    if (!(equal (exp_answer, real_answer))) {
    //if (!(m_traits.equal_2_object() (exp_answer, real_answer))) {
      if (!end_of_line_printed) print_end_of_line();
      std::cout << "Expected " << str << ": " << exp_answer << std::endl
                << "Obtained " << str << ": " << real_answer << std::endl;
      end_of_line_printed = true;
      return false;
    }
    return true;
  }

  bool compare(const unsigned int & exp_answer,
                    const unsigned int & real_answer,
                    const char * str = "result")
  {
    if (exp_answer != real_answer) {
      if (!end_of_line_printed) print_end_of_line();
      std::cout << "Expected " << str << ": " << exp_answer << std::endl
                << "Obtained " << str << ": " << real_answer << std::endl;
      end_of_line_printed = true;
      return false;
    }
    return true;
  }

  bool compare_and_print_for_points(const Point_2 & exp_answer, const Point_2 & real_answer,
                         const char * str = "point")
  {
    bool result = compare_points(exp_answer, real_answer, str);
    print_result(result);
    return result;
  }

  bool compare_and_print(const unsigned int & exp_answer, const unsigned int & real_answer,
                         const char * str = "result")
  {
    bool result = compare(exp_answer, real_answer, str);
    print_result(result);
    return result;
  }

  template <class stream>
  bool read_point(stream & is, Point_2 &);

  template <class stream>
  bool read_xcurve(stream & is, X_monotone_curve_2 &);

  template <class stream>
  bool read_curve(stream & is, Curve_2 &);

  //@{

  /*! Test Compare_x_2
   */
  bool compare_x_wrapper(std::istringstream & );
  bool compare_x_wrapper_imp(std::istringstream & , CGAL::Tag_false);
  bool compare_x_wrapper_imp(std::istringstream & , CGAL::Tag_true);

  /*! Test Boundary_in_x_2
   */
  bool boundary_in_x_wrapper(std::istringstream & );
  bool boundary_in_x_wrapper_imp(std::istringstream & , CGAL::Tag_false);
  bool boundary_in_x_wrapper_imp(std::istringstream & , CGAL::Tag_true);

  /*! Test Boundary_in_y_2
   */
  bool boundary_in_y_wrapper(std::istringstream & );
  bool boundary_in_y_wrapper_imp(std::istringstream & , CGAL::Tag_false);
  bool boundary_in_y_wrapper_imp(std::istringstream & , CGAL::Tag_true);

  /*! Compare_xy_2
   */
  bool compare_xy_wrapper(std::istringstream & line);

  /*! Tests Construct_min_vertex_2.
   *  Degenerate case: vertical curve.
   */
  bool min_vertex_wrapper(std::istringstream & line);

  /*! Tests Construct_max_vertex_2.
   * Degenerate case: vertical curve.
   */
  bool max_vertex_wrapper(std::istringstream & line);

  bool is_vertical_wrapper(std::istringstream & line);

  /*! Tests Compare_y_at_x_2.
   * Return the location of the given point with respect to the input curve.
   * Degenerate cases: The point is an endpoint of the curve.
   *                   The curve is vertical.
   */
  bool compare_y_at_x_wrapper(std::istringstream & line);

  /*! Tests Compare_y_at_x_left_2.
   * Compare the y value of two x-monotone curves immediately to the left
   * (resp. right) of their intersection point.
   * Degenerate cases: The curves coincide.
   *                   The curves coincide and vertical.
   *                   One of the curves is vertical.
   */
  bool compare_y_at_x_left_wrapper(std::istringstream & line);
  bool compare_y_at_x_left_wrapper_imp(std::istringstream & line,
                                       CGAL::Tag_false);
  bool compare_y_at_x_left_wrapper_imp(std::istringstream & line,
                                       CGAL::Tag_true);

  /*! Tests Compare_y_at_x_right_2.
   * Compare the y value of two x-monotone curves immediately to the right
   * (resp. right) of their intersection point.
   * Degenerate cases: The curves coincide.
   *                   The curves coincide and vertical.
   *                   One of the curves is vertical.
   */
  bool compare_y_at_x_right_wrapper(std::istringstream & line);

  /*! Tests Equal_2::operator()(Point_2, Point_2).
   * Check whether two points are the same.
   */
  bool equal_points_wrapper(std::istringstream & line);

  /*! Tests Equal_2::operator()(X_monotone_curve_2, X_monotone_curve_2).
   * Check whether two x-monotone curves are the same.
   */
  bool equal_curves_wrapper(std::istringstream & line);

  /*! Tests Make_x_monotone_2.
   * Cut the given curve into x-monotone subcurves and insert them into the
   * given output iterator.
   * Degenerate cases for polylines: The first segment is vertical. The last
   * segment is vertical. Both firt and last are vertical. An internal segment
   * is vertical.
   */
  bool make_x_monotone_wrapper(std::istringstream & line);

  /*! Tests Intersect_2.
   * Find the intersections of the two given curves and insert them into the 
   * given output iterator.
   * Degenerate cases for polylines: The most right (resp. left) endpoints of
   * the two curves coincide. Both endpoints coincide. The most right (resp.
   * left) endpoint of one curve and the first (resp. last) segment of the
   * other coincide.
   */
  bool intersect_wrapper(std::istringstream & line);

  /*! Tests Split_2.
   * Split a given x-monotone curve at a given point into two sub-curves.
   * Degenerate cases for polylines: the point and a polyline internal point
   * coincides.
   */
  bool split_wrapper(std::istringstream & line);

  /*! Tests Are_mergeable_2.
   * Check whether it is possible to merge two given x-monotone curves.
   */
  bool are_mergeable_wrapper(std::istringstream & line);
  bool are_mergeable_wrapper_imp(std::istringstream & line, CGAL::Tag_false);
  bool are_mergeable_wrapper_imp(std::istringstream & line, CGAL::Tag_true);

  /*! Tests Merge_2.
   * Merge two given x-monotone curves into a single curve.
   */
  bool merge_wrapper(std::istringstream & line);
  bool merge_wrapper_imp(std::istringstream & line, CGAL::Tag_false);
  bool merge_wrapper_imp(std::istringstream & line, CGAL::Tag_true);

  /*! Tests Approximate_2.
   * Return an approximation of a point coordinate.
   */
  bool approximate_wrapper(std::istringstream & line);

  /*! tests Construct_x_monotone_curve_2.
   * Return an x-monotone curve connecting the two given endpoints.
   */
  bool construct_x_monotone_curve_wrapper(std::istringstream & line);

  //@}

public:
  /*! Constructor */
  Traits_test(int argc, char * argv[]);

  /*! Destructor */
  ~Traits_test();

  /*! Entry point */
  bool start();
};

/*!
 * Constructor. 
 * Accepts test data file name.
 */
template <class T_Traits>
Traits_test<T_Traits>::Traits_test(int argc, char * argv[])
{
  typedef T_Traits Traits;
  abort_on_error=false;//run all tests (defualt)
  //abort_on_error=true;//run only one error test
  end_of_line_printed = true;
  violation_map[PRECONDITION]=std::string("precondition");
  violation_map[POSTCONDITION]=std::string("postcondition");
  violation_map[ASSERTION]=std::string("assertion");
  violation_map[WARNING]=std::string("warning");
  if (argc != 6) {
    std::cout << "Usage: " << argv[0] <<
      " points_file xcurves_file curves_file commands_file traits_type_name"
      << std::endl;
    end_of_line_printed = true;
  }
  else
  {
    m_filename_points = argv[1];
    m_filename_xcurves = argv[2];
    m_filename_curves = argv[3];
    m_filename_commands = argv[4];
    m_traitstype = argv[5];
  }
  m_wrappers[std::string("compare_x")] =
    &Traits_test<Traits>::compare_x_wrapper;
  m_wrappers[std::string("compare_xy")] =
    &Traits_test<Traits>::compare_xy_wrapper;
  m_wrappers[std::string("boundary_in_x")] =
    &Traits_test<Traits>::boundary_in_x_wrapper;
  m_wrappers[std::string("boundary_in_y")] =
    &Traits_test<Traits>::boundary_in_y_wrapper;
  m_wrappers[std::string("min_vertex")] =
    &Traits_test<Traits>::min_vertex_wrapper;
  m_wrappers[std::string("max_vertex")] =
    &Traits_test<Traits>::max_vertex_wrapper;
  m_wrappers[std::string("is_vertical")] =
    &Traits_test<Traits>::is_vertical_wrapper;
  m_wrappers[std::string("compare_y_at_x")] =
    &Traits_test<Traits>::compare_y_at_x_wrapper;
  m_wrappers[std::string("compare_y_at_x_left")] =
    &Traits_test<Traits>::compare_y_at_x_left_wrapper;
  m_wrappers[std::string("compare_y_at_x_right")] =
    &Traits_test<Traits>::compare_y_at_x_right_wrapper;
  m_wrappers[std::string("equal_points")] =
    &Traits_test<Traits>::equal_points_wrapper;
  m_wrappers[std::string("equal_curves")] =
    &Traits_test<Traits>::equal_curves_wrapper;
  m_wrappers[std::string("make_x_monotone")] =
    &Traits_test<Traits>::make_x_monotone_wrapper;
  m_wrappers[std::string("intersect")] =
    &Traits_test<Traits>::intersect_wrapper;
  m_wrappers[std::string("split")] =
    &Traits_test<Traits>::split_wrapper;
  m_wrappers[std::string("are_mergeable")] =
    &Traits_test<Traits>::are_mergeable_wrapper;
  m_wrappers[std::string("merge")] =
    &Traits_test<Traits>::merge_wrapper;
  m_wrappers[std::string("approximate")] =
    &Traits_test<Traits>::approximate_wrapper;
  m_wrappers[std::string("construct_x_monotone_curve")] =
    &Traits_test<Traits>::construct_x_monotone_curve_wrapper;
}

/*!
 * Destructor. 
 * Declares as virtual.
 */
template <class T_Traits>
Traits_test<T_Traits>::~Traits_test()
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
bool Traits_test<T_Traits>::start()
{
  std::ifstream in_pt(m_filename_points.c_str());
  std::ifstream in_xcv(m_filename_xcurves.c_str());
  std::ifstream in_cv(m_filename_curves.c_str());
  std::ifstream in_com(m_filename_commands.c_str());
  if (!in_pt.is_open()) {
    std::cerr << "Error opening file " << m_filename_points.c_str() << std::endl;
    end_of_line_printed = true;
    return false;
  }
  if (!in_xcv.is_open()) {
    std::cerr << "Error opening file " << m_filename_xcurves.c_str() << std::endl;
    end_of_line_printed = true;
    return false;
  }
  if (!in_cv.is_open()) {
    std::cerr << "Error opening file " << m_filename_curves.c_str() << std::endl;
    end_of_line_printed = true;
    return false;
  }
  if (!in_com.is_open()) {
    std::cerr << "Error opening file " << m_filename_commands.c_str() << std::endl;
    end_of_line_printed = true;
    return false;
  }
  if (!read_input(in_pt,POINT)) {
    in_pt.close(); 
    return false;
  }
  if (!read_input(in_xcv,XCURVE)) {
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
bool Traits_test<T_Traits>::read_input(std::ifstream & is ,Read_object_type obj_t)
{
  char one_line[128];
  skip_comments(is, one_line);
  std::istringstream str_stream(one_line, std::istringstream::in);
    try
    {
      for (int i = 0; !is.eof() ; ++i) 
      {
	  std::cout << " str_stream gcount = " << str_stream.gcount() << std::endl ;
          switch (obj_t)
          {
            case POINT :
	      m_points.resize(m_points.size()+1);
              if (!read_point(str_stream, m_points[i]))
              {
                 std::cerr << "Error reading point!" << std::endl;
                 end_of_line_printed = true;
                 return false;
               }
              break;
            case CURVE :
	      m_curves.resize(m_curves.size()+1);
              if (!read_curve(str_stream, m_curves[i]))
               {
                 std::cerr << "Error reading curves!" << std::endl;
                 end_of_line_printed = true;
                 return false;
              }
              break;
            case XCURVE :
	      m_xcurves.resize(m_xcurves.size()+1);
              if (!read_xcurve(str_stream, m_xcurves[i]))
              {
                 std::cerr << "Error reading xcurves!" << std::endl;
                 end_of_line_printed = true;
                 return false;
              }
              break;
	     }//switch
       str_stream.clear();
       skip_comments(is, one_line);
       str_stream.str(one_line);
     }//for
     /*for (int i = 0; i < num; ++i)
     { //for debug
        switch (obj_t)
        {
          case POINT :
            std::cout << "POINT size " << i << " " << m_points.size() << " := " << m_points[i] << std::endl;
            break;
          case CURVE :
            std::cout << "CURVE size " << i << " " << m_curves.size() << " := " << m_curves[i] << std::endl;
            break;
          case XCURVE :
            std::cout << "XCURVE size " << i << " " << m_xcurves.size() << " := " << m_xcurves[i] << std::endl;
            break;
        } //switch
     } */ //for
   }//try
   catch (std::exception e)
   {
      switch (obj_t)
      {
         case POINT :
           std::cout << "Abort Test, a violation during reading point!" 
                     << std::endl;
           break;
         case CURVE :
           std::cout << "Abort Test, a violation during reading curve!" 
                     << std::endl;
           break;
         case XCURVE :
           std::cout << "Abort Test, a violation during reading xcurve!" 
                     << std::endl;
           break;
      }//switch
      is.close();
      return false;
    }//catch
  return true;
}

/*!
 * Command dispatcher. Retrieves a line from the input file and performes
 * some action. See comments for suitable function in order to know specific
 * command arguments. 
 */
template <class T_Traits>
bool Traits_test<T_Traits>::perform(std::ifstream & is)
{
  bool test_result = true;
  std::cout << "Performing test : traits type is " << m_traitstype <<
    ", input files are " << m_filename_points << " " <<
    m_filename_xcurves << " " << m_filename_curves << " " <<
    m_filename_commands << std::endl;
  end_of_line_printed = true;
  char one_line[128];
  char buff[128];
  bool abort=false;
  int counter=0;
  while (!(is.eof() || abort)) 
  {
    skip_comments(is, one_line);
    std::istringstream str_stream(one_line, std::istringstream::in);
    buff[0] = '\0';
    str_stream.getline(buff, 128, ' ');
    std::string str_command(buff);
    unsigned int location=0;
    violation_occurred=violation_tested=NON;
    if ((int)str_command.find("_precondition",0)!=-1)
    {
      location=str_command.find("_precondition",0);
      violation_tested=PRECONDITION;
    }
    else if ((int)str_command.find("_postcondition",0)!=-1)
    {
      location=str_command.find("_postcondition",0);
      violation_tested=POSTCONDITION;
    }
    else if ((int)str_command.find("_assertion",0)!=-1)
    {
      location=str_command.find("_assertion",0);
      violation_tested=ASSERTION;
    }
    else if ((int)str_command.find("_warning",0)!=-1)
    {
      location=str_command.find("_warning",0);
      violation_tested=WARNING;
    }
    if (violation_tested!=NON)
    {
      str_command=str_command.substr(0,location);
      std::cout << "Test " << violation_map[violation_tested] 
                << " violation : ";
    }
    /*if (!test_result)
    std::cout << "bug" << std::endl;*/
    counter++;
    std::cout << "iter number : " << counter << std::endl;
    Wrapper_iter wi = m_wrappers.find(str_command);
    str_stream.clear();
    if (wi == m_wrappers.end()) continue;
    Wrapper wrapper = (*wi).second;
    try
    {
      test_result &= (this->*wrapper)(str_stream);
    }
    catch (Test_exception e)
    {
      Exception_type e_t = e.get();
      if ( e_t==UNEXPECTED_CONTINUE || e_t==UNEXPECTED_ABORT )
      {
        test_result=false;
      }
      else
      {
         //true for more information, and false for less
         bool display_all_violation_info=true;
         if (display_all_violation_info)
         {
            std::cout << "library " << e.library() << std::endl;
            std::cout << "filename " << e.filename() << std::endl;
            std::cout << "line_number " << e.line_number() << std::endl;
            std::cout << "expression " << e.expression() << std::endl;
         }
      }
      switch ( e_t )
      {
        case EXPECTED_CONTINUE :
          break;
        case EXPECTED_ABORT :
          std::cout << "Test successfully aborted" <<
            ", input files are " << m_filename_points << " " <<
            m_filename_xcurves << " " << m_filename_curves << " " <<
            m_filename_commands << std::endl;
          abort=true;
          break;
        case UNEXPECTED_CONTINUE :
          if (violation_occurred!=NON)
          {
            std::cout << "Unexpected " << violation_map[violation_occurred] <<
                         " violation occurred! next test ... " << std::endl; 
          }
          else
          {
            std::cout << "Violation did not occur!, expected " << 
                         violation_map[violation_tested] << 
                         " next test ... " << std::endl; 
          }
          break;
        case UNEXPECTED_ABORT :
          if (violation_occurred!=NON)
          {
            std::cout << "Unexpected " << violation_map[violation_occurred] <<
                         " violation occurred! abort ... "; 
          }
          else
          {
            std::cout << "Violation did not occur!, expected " << 
                         violation_map[violation_tested] << 
                         " abort ... "; 
          }
          std::cout << "input files are " << m_filename_points << " " <<
            m_filename_xcurves << " " << m_filename_curves << " " <<
            m_filename_commands << std::endl;
          abort=true;
          break;
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
void Traits_test<T_Traits>::skip_comments(std::ifstream & is, char * one_line)
{
  while (!is.eof()) {
    is.getline(one_line, 128);
    if (one_line[0] != '#') break;
  }
}

/*!
 */
template <class T_Traits>
std::string Traits_test<T_Traits>::remove_blanks(char * str)
{
  std::string result = "";
  bool flag = false;
  //only alphanumeric characters and underscores are allowed
  for (; *str != '\0'; ++str)
  {
    if ((*str >= '0' && *str <= '9') || //digits
        (*str >= 'A' && *str <= 'Z') || //upper case letters
        (*str >= 'a' && *str <= 'z') || //lower case letters
         *str=='_') //underscores
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
bool Traits_test<T_Traits_class>::translate_boolean(std::string & str_value)
{
  if (str_value == "TRUE") return true;
  return false;
}

/*!
 */
template <class T_Traits>
unsigned int
Traits_test<T_Traits>::translate_enumerator(std::string & str_value)
{
  if (str_value == "LARGER" ) {
    return static_cast<unsigned int>(CGAL::LARGER);
  } else if (str_value == "SMALLER" ) {
    return static_cast<unsigned int>(CGAL::SMALLER);
  } else if (str_value == "EQUAL" ) {
    return static_cast<unsigned int>(CGAL::EQUAL);
  }
  CGAL_assertion(false);
  return static_cast<unsigned int>(-220776); // My birthday :-)
}

/*!
 */
template <class T_Traits>
std::pair<Enum_type , unsigned int>
Traits_test<T_Traits>::translate_int_or_text(std::string & str_value)
{
  if (str_value == "MIN_END" ) {
    return std::pair<enum Enum_type,unsigned int>(CURVE_END,CGAL::MIN_END);
  } else if (str_value == "MAX_END" ) {
    return std::pair<enum Enum_type,unsigned int>(CURVE_END,CGAL::MAX_END);
  } else if (str_value == "SMALLER" ) {
    return std::pair<enum Enum_type,unsigned int>(SIGN,
                                    static_cast<unsigned int>(CGAL::SMALLER));
  } else if (str_value == "EQUAL" ) {
    return std::pair<enum Enum_type,unsigned int>(SIGN,
                                    static_cast<unsigned int>(CGAL::EQUAL));
  } else if (str_value == "LARGER" ) {
    return std::pair<enum Enum_type,unsigned int>(SIGN,
                                    static_cast<unsigned int>(CGAL::LARGER));
  } else if (str_value == "MINUS_INFINITY" ) {
    return std::pair<enum Enum_type,unsigned int>
      (BOUNDARY,static_cast<unsigned int>(CGAL::MINUS_INFINITY));
  } else if (str_value == "NO_BOUNDARY" ) {
    return std::pair<enum Enum_type,unsigned int>
      (BOUNDARY,static_cast<unsigned int>(CGAL::NO_BOUNDARY));
  } else if (str_value == "PLUS_INFINITY" ) {
    return std::pair<enum Enum_type,unsigned int>
      (BOUNDARY,static_cast<unsigned int>(CGAL::PLUS_INFINITY));
  }
  return std::pair<enum Enum_type,unsigned int>
    (NUMBER,static_cast<unsigned int>(atoi(str_value.c_str())));
}

/*!
 */
template <class T_Traits>
bool
Traits_test<T_Traits>::get_expected_boolean(std::istringstream & str_stream)
{
  char buff[128];
  str_stream.getline( buff, 128, '.');
  buff[str_stream.gcount()] = '\0';
  std::string str_expres = remove_blanks(buff);
  return translate_boolean(str_expres);
}

/*!
 */
template <class T_Traits>
unsigned int
Traits_test<T_Traits>::get_expected_enum(std::istringstream & str_stream)
{
  char buff[128];
  str_stream.getline(buff, 128, '.');
  buff[str_stream.gcount()] = '\0';
  std::string str_expres = remove_blanks(buff);
  return translate_enumerator(str_expres);
}

/*!
 */
template <class T_Traits>
std::pair<Enum_type,unsigned int>
Traits_test<T_Traits>::get_next_input(std::istringstream & str_stream)
{
  char buff[128];
  do
  {
    str_stream.getline(buff, 128, ' ');
  } while (str_stream.gcount()==1);
  buff[str_stream.gcount()] = '\0';
  std::string str_expres = remove_blanks(buff);
  return translate_int_or_text(str_expres);
}

/*! Test Boundary_in_x_2
*/
template <class T_Traits>
bool Traits_test<T_Traits>::boundary_in_x_wrapper(std::istringstream & str_stream)
{
  typedef typename T_Traits::Has_boundary_category       Has_boundary_category;
  return boundary_in_x_wrapper_imp(str_stream, Has_boundary_category());
}

template <class T_Traits>
bool Traits_test<T_Traits>::boundary_in_x_wrapper_imp(std::istringstream & , CGAL::Tag_false)
{
  CGAL_assertion(false);
  return false;
}

template <class T_Traits>
bool Traits_test<T_Traits>::boundary_in_x_wrapper_imp(std::istringstream & str_stream, CGAL::Tag_true)
{
  unsigned int id;
  str_stream >> id;
  std::pair<Enum_type,unsigned int> next_input = get_next_input(str_stream);
  if (next_input.first==CURVE_END)
  {
    CGAL::Curve_end cv_end=static_cast<CGAL::Curve_end>(next_input.second);
    next_input = get_next_input(str_stream);
    if (next_input.first==BOUNDARY)
    {
      CGAL::Boundary_type exp_answer = static_cast<CGAL::Boundary_type>(next_input.second);
      CGAL::Boundary_type real_answer = m_traits.boundary_in_x_2_object()(m_xcurves[id],cv_end);
      did_violation_occur();
      std::cout << "Test: boundary_in_x( " << m_xcurves[id] << " , "
                << (cv_end==CGAL::MIN_END?"MIN_END":"MAX_END") << " ) ? "
                << (exp_answer==CGAL::MINUS_INFINITY?"MINUS_INFINITY":
                   (exp_answer==CGAL::PLUS_INFINITY?"PLUS_INFINITY":"NO_BOUNDARY")) << " ";
      return compare_and_print(exp_answer, real_answer);
    }
  }
  //should not get here
  CGAL_assertion(false);
  return false;
}

/*! Test Boundary_in_y_2
*/
template <class T_Traits>
bool Traits_test<T_Traits>::boundary_in_y_wrapper(std::istringstream & str_stream)
{
  typedef typename T_Traits::Has_boundary_category       Has_boundary_category;
  return boundary_in_y_wrapper_imp(str_stream, Has_boundary_category());
}

template <class T_Traits>
bool Traits_test<T_Traits>::boundary_in_y_wrapper_imp(std::istringstream &, CGAL::Tag_false)
{
  CGAL_assertion(false);
  return false;
}

template <class T_Traits>
bool Traits_test<T_Traits>::boundary_in_y_wrapper_imp(std::istringstream & str_stream, CGAL::Tag_true)
{
  unsigned int id;
  str_stream >> id;
  std::pair<Enum_type,unsigned int> next_input = get_next_input(str_stream);
  if (next_input.first==CURVE_END)
  {
    CGAL::Curve_end cv_end=static_cast<CGAL::Curve_end>(next_input.second);
    next_input = get_next_input(str_stream);
    if (next_input.first==BOUNDARY)
    {
      CGAL::Boundary_type exp_answer = static_cast<CGAL::Boundary_type>(next_input.second);
      CGAL::Boundary_type real_answer = m_traits.boundary_in_y_2_object()(m_xcurves[id],cv_end);
      did_violation_occur();
      std::cout << "Test: boundary_in_x( " << m_xcurves[id] << " , "
                << (cv_end==CGAL::MIN_END?"MIN_END":"MAX_END") << " ) ? "
                << (exp_answer==CGAL::MINUS_INFINITY?"MINUS_INFINITY":
                   (exp_answer==CGAL::PLUS_INFINITY?"PLUS_INFINITY":"NO_BOUNDARY")) << " ";
      return compare_and_print(exp_answer, real_answer);
    }
  }
  //should not get here
  CGAL_assertion(false);
  return false;
}

/*! Test Compare_x_2
 */
template <class T_Traits>
bool Traits_test<T_Traits>::compare_x_wrapper(std::istringstream & str_stream)
{
  typedef typename T_Traits::Has_boundary_category       Has_boundary_category;
  return compare_x_wrapper_imp(str_stream, Has_boundary_category());
}

template <class T_Traits>
bool Traits_test<T_Traits>::compare_x_wrapper_imp
                            (std::istringstream & str_stream, CGAL::Tag_true)
{
  unsigned int id1, id2;
  unsigned int real_answer = 0, exp_answer = 0;
  str_stream >> id1;
  std::pair<Enum_type,unsigned int> exp_answer_1 = get_next_input(str_stream);
  if (exp_answer_1.first==NUMBER)
  {
    id2=exp_answer_1.second;
    std::pair<Enum_type,unsigned int> exp_answer_2 =get_next_input(str_stream);
    if (exp_answer_2.first==SIGN)
    {
      exp_answer = exp_answer_2.second;
      std::cout << "Test: compare_x( " << m_points[id1]
                << ", " << m_points[id2] << " ) ? " << exp_answer << " ";
      real_answer =m_traits.compare_x_2_object()(m_points[id1], m_points[id2]);
    }
    if (exp_answer_2.first==CURVE_END)
    {
      CGAL::Curve_end cv_end_1 = static_cast<CGAL::Curve_end>
                                            (exp_answer_2.second);
      exp_answer = get_expected_enum(str_stream);
      std::cout << "Test: compare_x( " << m_points[id1]
                << ", " << m_xcurves[id2] << ", "
                << (exp_answer_2.second==CGAL::MIN_END?" MIN_END ":" MAX_END ")
                << " ) ? " << exp_answer << " ";
        real_answer = m_traits.compare_x_2_object()
                               (m_points[id1], m_xcurves[id2],cv_end_1);
    }
  }
  else if (exp_answer_1.first==CURVE_END)
  {
    CGAL::Curve_end cv_end_1=static_cast<CGAL::Curve_end>(exp_answer_1.second);
    str_stream >> id2;
    std::pair<Enum_type,unsigned int> exp_answer_2 =get_next_input(str_stream);
    CGAL_assertion(exp_answer_2.first==CURVE_END);
    CGAL::Curve_end cv_end_2=static_cast<CGAL::Curve_end>(exp_answer_2.second);
    exp_answer = get_expected_enum(str_stream);
    std::cout << "Test: compare_x( " << m_xcurves[id1] << ", "
              << (exp_answer_1.second==CGAL::MIN_END ?" MIN_END ":" MAX_END ")
              << ", " << m_xcurves[id2] << ", "
              << (exp_answer_2.second==CGAL::MIN_END ?" MIN_END ":" MAX_END ")
              << " ) ? " << exp_answer << " ";
      real_answer = m_traits.compare_x_2_object()(m_xcurves[id1], cv_end_1,
                                                  m_xcurves[id2], cv_end_2);
  }
  did_violation_occur();
  return compare_and_print(exp_answer, real_answer);
}

template <class T_Traits>
bool Traits_test<T_Traits>::compare_x_wrapper_imp
                            (std::istringstream & str_stream, CGAL::Tag_false)
{
  unsigned int id1, id2;
  str_stream >> id1 >> id2;
  unsigned int exp_answer = get_expected_enum(str_stream);
  unsigned int real_answer =
  m_traits.compare_x_2_object()(m_points[id1], m_points[id2]);
  did_violation_occur();
  std::cout << "Test: compare_x( " << m_points[id1] << ", "
            << m_points[id2] << " ) ? " << exp_answer << " ";
 return compare_and_print(exp_answer, real_answer);
}

/*! Test Compare_xy_2
 */
template <class T_Traits>
bool Traits_test<T_Traits>::compare_xy_wrapper(std::istringstream & str_stream)
{
  unsigned int id1, id2;
  str_stream >> id1 >> id2;
  unsigned int exp_answer = get_expected_enum(str_stream);
  unsigned int real_answer =
    m_traits.compare_xy_2_object()(m_points[id1], m_points[id2]);
  did_violation_occur();
  std::cout << "Test: compare_xy( " << m_points[id1] << ", "
            << m_points[id2] << " ) ? ";
  std::cout << exp_answer << " ";
  return compare_and_print(exp_answer, real_answer);
}

/*! Tests Construct_min_vertex_2.
 * Degenerate case: vertical curve.
 */
template <class T_Traits>
bool Traits_test<T_Traits>::min_vertex_wrapper(std::istringstream & str_stream)
{
  unsigned int id1, id2;
  str_stream >> id1 >> id2;
  Point_2 & exp_answer = m_points[id2];
  Point_2 real_answer =
    m_traits.construct_min_vertex_2_object()(m_xcurves[id1]);
  did_violation_occur();
  std::cout << "Test: min_vertex( " << m_xcurves[id1] << " ) ? ";
  std::cout << exp_answer << " ";
  return compare_and_print_for_points(exp_answer, real_answer);
}

/*! Tests Construct_max_vertex_2.
 * Degenerate case: vertical curve.
 */
template <class T_Traits>
bool Traits_test<T_Traits>::max_vertex_wrapper(std::istringstream & str_stream)
{
  unsigned int id1, id2;
  str_stream >> id1 >> id2;
  Point_2 & exp_answer = m_points[id2];
  Point_2 real_answer =
    m_traits.construct_max_vertex_2_object()(m_xcurves[id1]);
  did_violation_occur();
  std::cout << "Test: max_vertex( " << m_xcurves[id1] << " ) ? ";
  std::cout << exp_answer << " ";
  return compare_and_print_for_points(exp_answer, real_answer);
}

template <class T_Traits>
bool
Traits_test<T_Traits>::is_vertical_wrapper(std::istringstream & str_stream)
{
  unsigned int id;
  str_stream >> id;
  bool exp_answer = get_expected_boolean(str_stream);
  bool real_answer = m_traits.is_vertical_2_object()(m_xcurves[id]);
/*  std::cout << "id " << id << (exp_answer?" exp_answer true ":" exp_answer false ")
            << (real_answer?" real_answer true ":" real_answer false ") << std::endl;*/
  did_violation_occur();
  std::cout << "Test: is_vertical( " << m_xcurves[id]
            << " ) ? ";
  std::cout << exp_answer << " ";
  return compare_and_print(exp_answer, real_answer);
}

/*! Tests Compare_y_at_x_2.
 * Return the location of the given point with respect to the input curve.
 * Degenerate cases: The point is an endpoint of the curve.
 *                   The curve is vertical.
 */
template <class T_Traits>
bool
Traits_test<T_Traits>::compare_y_at_x_wrapper(std::istringstream & str_stream)
{
  unsigned int id1, id2;
  str_stream >> id1 >> id2;
  unsigned int exp_answer = get_expected_enum(str_stream);
  unsigned int real_answer =
    m_traits.compare_y_at_x_2_object()(m_points[id1], m_xcurves[id2]);
  did_violation_occur();
  std::cout << "Test: compare_y_at_x( " << m_points[id1] << ","
            << m_xcurves[id2]
            << " ) ? ";
  std::cout << exp_answer << " ";
  return compare_and_print(exp_answer, real_answer);
}

/*! Tests Compare_y_at_x_left_2.
 * Compare the y value of two x-monotone curves immediately to the left
 * of their intersection point.
 * Degenerate cases: The curves coincide.
 *                   The curves coincide and vertical.
 *                   One of the curves is vertical.
 */
template <class T_Traits>
bool
Traits_test<T_Traits>::
compare_y_at_x_left_wrapper(std::istringstream & str_stream)
{
  typedef typename T_Traits::Has_left_category          Has_left_category;
  return compare_y_at_x_left_wrapper_imp(str_stream, Has_left_category());
}

template <class T_Traits>
bool
Traits_test<T_Traits>::
compare_y_at_x_left_wrapper_imp(std::istringstream & ,
                                CGAL::Tag_false)
{
//  std::istringstream dummy_stream(str_stream); //to avoid warnings of unused variable
  CGAL_assertion(false);
  return false;
}

template <class T_Traits>
bool
Traits_test<T_Traits>::
compare_y_at_x_left_wrapper_imp(std::istringstream & str_stream,
                                CGAL::Tag_true)
{
  unsigned int id1, id2, id3;
  str_stream >> id1 >> id2 >> id3;
  unsigned int exp_answer = get_expected_enum(str_stream);
  unsigned int real_answer = m_traits.compare_y_at_x_left_2_object()
                               (m_xcurves[id1],m_xcurves[id2],m_points[id3]);
  did_violation_occur();
  std::cout << "Test: compare_y_at_x_left( " << m_xcurves[id1] << ","
            << m_xcurves[id2] << ", " << m_points[id3]
            << " ) ? "; 
  std::cout << exp_answer << " ";

  return compare_and_print(exp_answer, real_answer);
}

/*! Tests Compare_y_at_x_right_2.
 * Compare the y value of two x-monotone curves immediately to the right
 * of their intersection point.
 * Degenerate cases: The curves coincide.
 *                   The curves coincide and vertical.
 *                   One of the curves is vertical.
 */
template <class T_Traits>
bool
Traits_test<T_Traits>::
compare_y_at_x_right_wrapper(std::istringstream & str_stream)
{
  unsigned int id1, id2, id3;
  str_stream >> id1 >> id2 >> id3;
  unsigned int exp_answer = get_expected_enum(str_stream);
  unsigned int real_answer =
    m_traits.compare_y_at_x_right_2_object()
      (m_xcurves[id1], m_xcurves[id2], m_points[id3]);
  did_violation_occur();
  std::cout << "Test: compare_y_at_x_right( " << m_xcurves[id1] << ","
            << m_xcurves[id2] << ", " << m_points[id3]
            << " ) ? ";
  std::cout << exp_answer << " ";
  return compare_and_print(exp_answer, real_answer);
}

/*! Tests Equal_2::operator()(Point_2, Point_2).
 * Check whether two points are the same.
 */
template <class T_Traits>
bool
Traits_test<T_Traits>::equal_points_wrapper(std::istringstream & str_stream)
{
  unsigned int id1, id2;
  str_stream >> id1 >> id2;
  bool exp_answer = get_expected_boolean(str_stream);
  bool real_answer = m_traits.equal_2_object()(m_points[id1], m_points[id2]);
  did_violation_occur();
  std::cout << "Test: equal( " << m_points[id1] << ", "
            << m_points[id2] << " ) ? ";
  std::cout << exp_answer << " ";
  return compare_and_print(exp_answer, real_answer);
}

/*! Tests Equal_2::operator()(X_monotone_curve_2, X_monotone_curve_2).
 * Check whether two x-monotone curves are the same.
 */
template <class T_Traits>
bool
Traits_test<T_Traits>::equal_curves_wrapper(std::istringstream & str_stream)
{
  unsigned int id1, id2;
  str_stream >> id1 >> id2;
  bool exp_answer = get_expected_boolean(str_stream);
  bool real_answer = m_traits.equal_2_object()(m_xcurves[id1], m_xcurves[id2]);
  did_violation_occur();
  std::cout << "Test: equal( " << m_xcurves[id1] << ", "
            << m_xcurves[id2] << " ) ? ";
  std::cout << exp_answer << " ";
  return compare_and_print(exp_answer, real_answer);
}

/*! Tests Make_x_monotone_2.
 * Cut the given curve into x-monotone subcurves and insert them into the
 * given output iterator.
 * Degenerate cases for polylines: The first segment is vertical. The last
 * segment is vertical. Both firt and last are vertical. An internal segment
 * is vertical.
 */
template <class T_Traits>
bool
Traits_test<T_Traits>::make_x_monotone_wrapper(std::istringstream & str_stream)
{
  typedef T_Traits                              Traits;
  typedef typename Traits::Point_2              Point_2;
  typedef typename Traits::X_monotone_curve_2   X_monotone_curve_2;
  typedef typename Traits::Curve_2              Curve_2;
  typedef typename Traits::Equal_2              Equal_2;
  
  unsigned int id;
  str_stream >> id;
  std::vector<CGAL::Object> object_vec;
  m_traits.make_x_monotone_2_object()(m_curves[id],
                                      std::back_inserter(object_vec));
  did_violation_occur();
  std::cout << "Test: make_x_monotone( " << m_curves[id]
            << " ) ? ";
  unsigned int num;
  str_stream >> num;
  if (!compare(num, (unsigned int)object_vec.size(), "size")) {
    print_result(false);
    return false;
  }
  Equal_2 equal = m_traits.equal_2_object();
  for (unsigned int i = 0; i < num; ++i) {
    unsigned int type;                  // 0 - point, 1 - x-monotone curve
    str_stream >> type;

    unsigned int id;                    // The id of the point or x-monotone
    str_stream >> id;                   // ... curve respectively

    unsigned int exp_type = 1;
    const X_monotone_curve_2 * xcv_ptr;
    xcv_ptr = CGAL::object_cast<X_monotone_curve_2> (&(object_vec[i]));
    if (xcv_ptr != NULL) {
      if (!compare(type, exp_type, "type")) {
        print_result(false);
        return false;
      }
      if (!equal(m_xcurves[id], *xcv_ptr)) {
        std::cerr << "Expected x-monotone curve [" << id << "]: "
                  << m_xcurves[id] << std::endl
                  << "Obtained x-monotone curve: "
                  << *xcv_ptr << std::endl;
        end_of_line_printed = true;
        print_result(false);
        return false;
      }
      continue;
    }

    exp_type = 0;
    const Point_2 * pt_ptr;
    pt_ptr = CGAL::object_cast<Point_2> (&(object_vec[i]));
    CGAL_assertion (pt_ptr != NULL);
    if (!compare(type, exp_type, "type")) {
      print_result(false);
      return false;
    }
    if (!equal(m_points[id], *pt_ptr)) {
      std::cerr << "Expected point [" << i << "]: "
                << m_points[id] << std::endl
                << "Obtained point: "
                << *pt_ptr << std::endl;
      end_of_line_printed = true;
      print_result(false);
      return false;
    }
  }
  object_vec.clear();
  print_result(true);
  return true;
}

/*! Tests Intersect_2.
 * Find the intersections of the two given curves and insert them into the 
 * given output iterator.
 * Degenerate cases for polylines: The most right (resp. left) endpoints of
 * the two curves coincide. Both endpoints coincide. The most right (resp.
 * left) endpoint of one curve and the first (resp. last) segment of the
 * other coincide.
 */
template <class T_Traits>
bool Traits_test<T_Traits>::intersect_wrapper(std::istringstream & str_stream)
{
  typedef T_Traits                              Traits;
  typedef typename Traits::Point_2              Point_2;
  typedef typename Traits::X_monotone_curve_2   X_monotone_curve_2;
  typedef typename Traits::Curve_2              Curve_2;
  typedef typename Traits::Equal_2              Equal_2;

  unsigned int id1, id2;
  str_stream >> id1 >> id2;
  std::vector<CGAL::Object> object_vec;
  (void) m_traits.intersect_2_object()(m_xcurves[id1], m_xcurves[id2],
                                std::back_inserter(object_vec));
  did_violation_occur();
  std::cout << "Test: intersect( " << m_xcurves[id1] << "," << m_xcurves[id2]
            << " ) ? ";
  unsigned int num;
  str_stream >> num;
  if (!compare(num, object_vec.size(), "size")) {
    print_result(false);
    return false;
  }
  Equal_2 equal = m_traits.equal_2_object();
  for (unsigned int i = 0; i < num; ++i) {
    unsigned int type;                  // 0 - point, 1 - x-monotone curve
    str_stream >> type;
    unsigned int id;                    // The id of the point or x-monotone
    str_stream >> id;                   // ... curve respectively
    unsigned int monotonicity;
    if (type == 0) {
      str_stream >> monotonicity;
    }
    unsigned int exp_type = 1;
    const X_monotone_curve_2 * xcv_ptr;
    xcv_ptr = CGAL::object_cast<X_monotone_curve_2> (&(object_vec[i]));
    if (xcv_ptr != NULL) {
      if (!compare(type, exp_type, "type")) {
        print_result(false);
        return false;
      }
      if (!equal(m_xcurves[id], *xcv_ptr)) {
        std::cerr << "Expected x-monotone curve [" << i << "]: "
                  << m_xcurves[id] << std::endl
                  << "Obtained x-monotone curve: "
                  << *xcv_ptr << std::endl;
        end_of_line_printed = true;
        print_result(false);
        return false;
      }
      continue;
    }

    exp_type = 0;
    typedef std::pair<Point_2,unsigned int> Point_2_pair;
    const Point_2_pair * pt_pair_ptr;
    pt_pair_ptr = CGAL::object_cast<Point_2_pair> (&(object_vec[i]));
    CGAL_assertion (pt_pair_ptr != NULL);
    if (!compare(type, exp_type, "type")) {
      print_result(false);
      return false;
    }
    if (!equal(m_points[id], (*pt_pair_ptr).first)) {
      std::cerr << "Expected point [" << i << "]: "
                << m_points[id] << std::endl
                << "Obtained point: "
                << (*pt_pair_ptr).first << std::endl;
      end_of_line_printed = true;
      print_result(false);
      return false;
    }
    if (!compare(monotonicity, (*pt_pair_ptr).second, "monotonicity")) {
      print_result(false);
      return false;
    }
  }
  object_vec.clear();
  print_result(true);
  return true;
}

/*! Tests Split_2.
 * Split a given x-monotone curve at a given point into two sub-curves.
 * Degenerate cases for polylines: the point and a polyline internal point
 * coincides.
 */
template <class T_Traits>
bool Traits_test<T_Traits>::split_wrapper(std::istringstream & str_stream)
{
  typedef T_Traits                              Traits;
  typedef typename Traits::Point_2              Point_2;
  typedef typename Traits::X_monotone_curve_2   X_monotone_curve_2;
  typedef typename Traits::Curve_2              Curve_2;
  typedef typename Traits::Equal_2              Equal_2;

  unsigned int id1, id2;
  str_stream >> id1 >> id2;
  X_monotone_curve_2 cv1, cv2;
  m_traits.split_2_object()(m_xcurves[id1], m_points[id2], cv1, cv2);
  did_violation_occur();
  std::cout << "Test: split( " << m_xcurves[id1] << "," << m_points[id2]
            << " ) ? ";
  Equal_2 equal = m_traits.equal_2_object();
  str_stream >> id1 >> id2;                   // ... curve respectively
  if (!equal(m_xcurves[id1], cv1) || !equal(m_xcurves[id2], cv2)) {
    std::cerr << "Expected x-monotone curves: "
              << m_xcurves[id1] << ", " << m_xcurves[id2] << std::endl
              << "Obtained x-monotone curve: "
              << cv1 << "," << cv2 << std::endl;
    end_of_line_printed = true;
    print_result(false);
    return false;
  }
  print_result(true);
  return true;
}

/*! Tests Are_mergeable_2.
 * Check whether it is possible to merge two given x-monotone curves.
 */
template <class T_Traits>
bool
Traits_test<T_Traits>::are_mergeable_wrapper(std::istringstream & str_stream)
{
  typedef typename T_Traits::Has_merge_category         Has_merge_category;
  return are_mergeable_wrapper_imp(str_stream, Has_merge_category());
}

template <class T_Traits>
bool
Traits_test<T_Traits>::
are_mergeable_wrapper_imp(std::istringstream & , CGAL::Tag_false)
{
//  std::istringstream dummy_stream = str_stream; //to avoid warnings of unused variable
  CGAL_assertion(false);
  return false;
}

template <class T_Traits>
bool
Traits_test<T_Traits>::
are_mergeable_wrapper_imp(std::istringstream & str_stream, CGAL::Tag_true)
{
  unsigned int id1, id2;
  str_stream >> id1 >> id2;
  bool exp_answer = get_expected_boolean(str_stream);
  bool real_answer = m_traits.are_mergeable_2_object()(m_xcurves[id1],
                                                       m_xcurves[id2]);
  did_violation_occur();
  std::cout << "Test: are_mergeable( " << m_xcurves[id1] << ", "
            << m_xcurves[id2] << " ) ? ";
  std::cout << exp_answer << " ";
  return compare_and_print(exp_answer, real_answer);
}

/*! Tests Merge_2.
 * Merge two given x-monotone curves into a single curve.
 */
template <class T_Traits>
bool Traits_test<T_Traits>::merge_wrapper(std::istringstream & str_stream)
{
  typedef typename T_Traits::Has_merge_category         Has_merge_category;
  return merge_wrapper_imp(str_stream, Has_merge_category());
}

template <class T_Traits>
bool Traits_test<T_Traits>::merge_wrapper_imp(std::istringstream & ,
                                              CGAL::Tag_false)
{
//  std::istringstream dummy_stream(str_stream); //to avoid warnings of unused variable
  CGAL_assertion(false);
  return false;
}

template <class T_Traits>
bool Traits_test<T_Traits>::merge_wrapper_imp(std::istringstream & str_stream,
                                              CGAL::Tag_true)
{
  typedef T_Traits                              Traits;
  typedef typename Traits::X_monotone_curve_2   X_monotone_curve_2;
  typedef typename Traits::Equal_2              Equal_2;

  unsigned int id1, id2, id;
  str_stream >> id1 >> id2 >> id;
  X_monotone_curve_2 cv;
  m_traits.merge_2_object()(m_xcurves[id1], m_xcurves[id2], cv);
  did_violation_occur();
  std::cout << "Test: merge( " << m_xcurves[id1] << ", "
            << m_xcurves[id2] << " ) ? " << m_xcurves[id] << " ";
  Equal_2 equal = m_traits.equal_2_object();
  if (!equal(m_xcurves[id], cv)) {
    std::cerr << "Expected x-monotone curve: " << m_xcurves[id] << std::endl
              << "Obtained x-monotone curve: " << cv << std::endl;
    end_of_line_printed = true;
    print_result(false);
    return false;
  }
  print_result(true);
  return true;
}

/*! Tests Approximate_2.
 * Return an approximation of a point coordinate.
 */
template <class T_Traits>
bool
Traits_test<T_Traits>::approximate_wrapper(std::istringstream & )
{
//  std::istringstream dummy_stream(str_stream); //to avoid warnings of unused variable
  return false;
}

/*! tests Construct_x_monotone_curve_2.
 * Return an x-monotone curve connecting the two given endpoints.
 */
template <class T_Traits>
bool Traits_test<T_Traits>::
construct_x_monotone_curve_wrapper(std::istringstream & )
{
//  std::istringstream dummy_stream(str_stream); //to avoid warnings of unused variable
  return false;
}

template <class T_Traits>
template <class stream>
bool
Traits_test<T_Traits>::read_point(stream & is, typename T_Traits::Point_2 & p)
{
  Basic_number_type x, y;
  is >> x >> y;
  p = typename T_Traits::Point_2(x, y);
  return true;
}

template <class T_Traits>
template <class stream>
bool
Traits_test<T_Traits>::read_xcurve(stream & is,
                                   typename T_Traits::X_monotone_curve_2 & xcv)
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
Traits_test<T_Traits>::read_curve(stream & is, typename T_Traits::Curve_2 & cv)
{
  Basic_number_type x1, y1, x2, y2;
  is >> x1 >> y1 >> x2 >> y2;
  Point_2 p1(x1, y1);
  Point_2 p2(x2, y2);
  CGAL_assertion(p1 != p2);
  cv = typename T_Traits::Curve_2(p1, p2);
  return true;
}

#endif
