#ifndef CGAL_TRAITS_TEST_H
#define CGAL_TRAITS_TEST_H

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>

template <class T_Traits>
class Traits_test {
private:
  typedef T_Traits                                      Traits;
  typedef typename Traits::Point_2                      Point_2;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename Traits::Curve_2                      Curve_2;

  /*! The input data file */
  std::string m_filename;

  /*! An instance of the traits */
  Traits m_traits;

  /*! The container of input points */
  std::vector<Point_2>  m_points;

  /*! The container of input curves */
  std::vector<Curve_2>  m_curves;

  /*! A map between (strings) commands and (member functions) operations */
  typedef bool (Traits_test::* Wrapper)(std::istringstream &);
  typedef std::map<std::string, Wrapper>        Wrapper_map;
  typedef typename Wrapper_map::iterator        Wrapper_iter;
  Wrapper_map m_wrappers;
  
  /*! Collect the data for the test */
  bool collect_data(std::ifstream & is);

  /*! Perform the test */
  bool perform(std::ifstream & is);
  
  /*! Skip comments */
  void skip_comments(std::ifstream & is, char * one_line);

  std::string remove_blanks(char * str);
  bool get_expected_boolean(std::istringstream & str_line);
  unsigned int get_expected_enum(std::istringstream & str_line);
  bool translate_boolean(std::string & str_value);
  unsigned int translate_enumerator(std::string & str_value);

  template <class T>  
  bool print_was_successful_or_not(T & exp_answer, T & real_answer)
  {
    bool ret_value = true;
    if (exp_answer == real_answer) {
      std::cout << "Passed"     << std::endl;
    } else{
      std::cout << "Failed" << std::endl
                << "Expected result: " << exp_answer << std::endl
                << "Obtained result: " << real_answer << std::endl;
      ret_value = false;
    } 
    return ret_value;      
  }

  bool read_curve(std::ifstream & is, Curve_2 & cv);
  bool read_point(std::ifstream & is, Point_2 & p);
  
  //@{

  /*! Test Compare_x_2
   */
  bool compare_x_wrapper(std::istringstream & line);

  /*! Compare_xy_2
   */
  bool compare_xy_wrapper(std::istringstream & line);  
  
  /*! Tests Construct_min_vertex_2.
   * Degenerate case: vertical curve.
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

  /*! Tests Compare_y_at_x_right_2.
   * Compare the y value of two x-monotone curves immediately to the right
   * (resp. right) of their intersection point.
   * Degenerate cases: The curves coincide.
   *                   The curves coincide and vertical.
   *                   One of the curves is vertical.
   */
  bool compare_y_at_x_right_wrapper(std::istringstream & line);
  
  /*! Tests Equal_2.
   * Check whether two points are the same.
   * Check whether two x-monotone curves are the same.
   */
  bool equal_wrapper(std::istringstream & line);

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

  /*! Tests Merge_2.
   * Merge two given x-monotone curves into a single curve.
   */
  bool merge_wrapper(std::istringstream & line);

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
  
  if (argc < 2 || argc >= 3)
    std::cout << "Usage: " << argv[0] << " test_data_file" << std::endl;  
  else
    m_filename = argv[1];

  m_wrappers[std::string("compare_x")] =
    &Traits_test<Traits>::compare_x_wrapper;
  m_wrappers[std::string("compare_xy")] =
    &Traits_test<Traits>::compare_xy_wrapper;  
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
  m_wrappers[std::string("equal")] =
    &Traits_test<Traits>::equal_wrapper;
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
  m_filename.clear();
  m_points.clear();
  m_curves.clear();
}

/*!
 * Test entry point 
 */
template<class T_Traits>
bool Traits_test<T_Traits>::start()
{
  std::ifstream is(m_filename.c_str());
  if (!is.is_open()) {
    std::cerr << "Error opening file " << m_filename.c_str();
    return false;
  }
  if (!collect_data(is)) {
    is.close(); 
    return false;
  }
  if (!perform(is)) {
    is.close(); 
    return false;
  }
  return true;
}

/*!
 * Collects data. Fills all_curves_vec and all_points_vec
 */
template <class T_Traits>
bool Traits_test<T_Traits>::collect_data(std::ifstream & is)
{
  char one_line[128];
  unsigned int n_curves, n_points;
  unsigned int i;

  skip_comments(is, one_line);
  std::string string_values1(one_line);
  std::istringstream str_line1(string_values1, std::istringstream::in);
  str_line1 >> n_curves;
  for (i = 0; i < n_curves; ++i) {
    Curve_2 cv;
    if (!read_curve(is, cv)) {
      std::cerr << "Error reading curve!" << std::endl;
      return false;
    }
    m_curves.push_back(cv);
  }
  skip_comments(is, one_line);
  std::string string_values2(one_line);
  std::istringstream str_line2(string_values2, std::istringstream::in);
  str_line2 >> n_points;
  for (i = 0; i < n_points; ++i) {
    Point_2 p;
    if (!read_point(is, p)) {
      std::cerr << "Error reading point!" << std::endl;
      return false;
    }
    m_points.push_back(p);
  }
  while (!is.eof() && is.get() != '\n');
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
  std::cout << "Performing test ..." << std::endl;  
  char one_line[128];
  char buff[128];
  while (!is.eof()) {
    skip_comments(is, one_line);
    buff[0] = '\0';
    std::string string_values(one_line);
    std::istringstream str_line(string_values, std::istringstream::in);
    str_line.getline(buff, 128, ' ');
    std::string str_command(buff);
    Wrapper_iter wi = m_wrappers.find(str_command);
    if (wi == m_wrappers.end()) continue;
    Wrapper wrapper = (*wi).second;
    test_result &= (this->*wrapper)(str_line);
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
    // std::cerr << std::endl; // SUNPRO likes this...      
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
  while (*str != '\0') {
    if (*str != ' ') result += *str;
    ++str;
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

  return static_cast<unsigned int>(-220776); // My birthday :-)
}

/*!
 */
template <class T_Traits>
bool Traits_test<T_Traits>::get_expected_boolean(std::istringstream & str_line)
{
  char buff[128];
  strLine.getline( buff, 128, '.' );
  buff[str_line.gcount()] = '\0';
  std::string str_expres = remove_blanks(buff);
  return translate_boolean(str_expres);
}

/*!
 */
template <class T_Traits>
unsigned int
Traits_test<T_Traits>::get_expected_enum(std::istringstream & str_line)
{
  char buff[128];
  str_line.getline(buff, 128, '.');
  buff[str_line.gcount()] = '\0';
  std::string str_expres = remove_blanks(buff);
  return translate_enumerator(str_expres);
}

/*! Test Compare_x_2
 */
template <class T_Traits>
bool Traits_test<T_Traits>::compare_x_wrapper(std::istringstream & str_line)
{
  unsigned int id1, id2;
  str_line >> id1 >> id2;
  unsigned int exp_answer = get_expected_enum(str_line);
  unsigned int real_answer =
    m_traits.compare_x_2_object()(m_points[id1], m_points[id2]);
  std::cout << "Test: compare_x( " << m_points[id1] << ", "
            << m_points[id2] << " ) ? " << exp_answer << std::endl;
  return print_was_successful_or_not(exp_answer, real_answer);
}

/*! Test Compare_xy_2
 */
template <class T_Traits>
bool Traits_test<T_Traits>::compare_xy_wrapper(std::istringstream & str_line)
{
  unsigned int id1, id2;
  str_line >> id1 >> id2;
  unsigned int exp_answer = get_expected_enum(str_line);
  unsigned int real_answer =
    m_traits.compare_xy_2_object()(m_points[id1], m_points[id2]);
  std::cout << "Test: compare_xy( " << m_points[id1] << ", "
            << m_points[id2] << " ) ? " << exp_answer << std::endl;
  return print_was_successful_or_not(exp_answer, real_answer);
}

/*! Tests Construct_min_vertex_2.
 * Degenerate case: vertical curve.
 */
template <class T_Traits>
bool Traits_test<T_Traits>::min_vertex_wrapper(std::istringstream & line)
{
  return false;
}

/*! Tests Construct_max_vertex_2.
 * Degenerate case: vertical curve.
 */
template <class T_Traits>
bool Traits_test<T_Traits>::max_vertex_wrapper(std::istringstream & line)
{
  return false;
}

template <class T_Traits>
bool Traits_test<T_Traits>::is_vertical_wrapper(std::istringstream & line)
{
  return false;
}

/*! Tests Compare_y_at_x_2.
 * Return the location of the given point with respect to the input curve.
 * Degenerate cases: The point is an endpoint of the curve.
 *                   The curve is vertical.
 */
template <class T_Traits>
bool Traits_test<T_Traits>::compare_y_at_x_wrapper(std::istringstream & line)
{
  return false;
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
Traits_test<T_Traits>::compare_y_at_x_left_wrapper(std::istringstream & line)
{
  return false;
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
Traits_test<T_Traits>::compare_y_at_x_right_wrapper(std::istringstream & line)
{
  return false;
}  

/*! Tests Equal_2.
 * Check whether two points are the same.
 * Check whether two x-monotone curves are the same.
 */
template <class T_Traits>
bool Traits_test<T_Traits>::equal_wrapper(std::istringstream & line)
{
  return false;
}

/*! Tests Make_x_monotone_2.
 * Cut the given curve into x-monotone subcurves and insert them into the
 * given output iterator.
 * Degenerate cases for polylines: The first segment is vertical. The last
 * segment is vertical. Both firt and last are vertical. An internal segment
 * is vertical.
 */
template <class T_Traits>
bool Traits_test<T_Traits>::make_x_monotone_wrapper(std::istringstream & line)
{
  return false;
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
bool Traits_test<T_Traits>::intersect_wrapper(std::istringstream & line)
{
  return false;
}

/*! Tests Split_2.
 * Split a given x-monotone curve at a given point into two sub-curves.
 * Degenerate cases for polylines: the point and a polyline internal point
 * coincides.
 */
template <class T_Traits>
bool Traits_test<T_Traits>::split_wrapper(std::istringstream & line)
{
  return false;
}

/*! Tests Are_mergeable_2.
 * Check whether it is possible to merge two given x-monotone curves.
 */
template <class T_Traits>
bool Traits_test<T_Traits>::are_mergeable_wrapper(std::istringstream & line)
{
  return false;
}

/*! Tests Merge_2.
 * Merge two given x-monotone curves into a single curve.
 */
template <class T_Traits>
bool Traits_test<T_Traits>::merge_wrapper(std::istringstream & line)
{
  return false;
}

/*! Tests Approximate_2.
 * Return an approximation of a point coordinate.
 */
template <class T_Traits>
bool Traits_test<T_Traits>::approximate_wrapper(std::istringstream & line)
{
  return false;
}

/*! tests Construct_x_monotone_curve_2.
 * Return an x-monotone curve connecting the two given endpoints.
 */
template <class T_Traits>
bool Traits_test<T_Traits>::
construct_x_monotone_curve_wrapper(std::istringstream & line)
{
  return false;
}

template <class T_Traits>
bool
Traits_test<T_Traits>::read_curve(std::ifstream & is,
                                  typename T_Traits::Curve_2 & cv)
{
  CGAL_assertion_msg(0, "Not specialized!!!");
  return false;
}

template <class T_Traits>
bool
Traits_test<T_Traits>::read_point(std::ifstream & is,
                                  typename T_Traits::Point_2 & p)
{
  CGAL_assertion_msg(0, "Not specialized!!!");
  return false;
}

#endif
