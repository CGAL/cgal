#ifndef CGAL_TRAITS_TEST_H
#define CGAL_TRAITS_TEST_H

#include <CGAL/Object.h>

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

  /*! The container of x-monotone curves */
  std::vector<X_monotone_curve_2>  m_xcurves;

  /*! Indicates whether the end-of-line has been printed */
  bool m_end_of_line_printed;
  
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
  bool get_expected_boolean(std::istringstream & str_stream);
  unsigned int get_expected_enum(std::istringstream & str_stream);
  bool translate_boolean(std::string & str_value);
  unsigned int translate_enumerator(std::string & str_value);

  void print_end_of_line()
  {
    std::cout << std::endl;
    m_end_of_line_printed = true;
  }
  
  void print_result(bool result)
  {
    if (!result && !m_end_of_line_printed) print_end_of_line();
    std::cout << ((result) ? "Passed" : "Failed");
    print_end_of_line();
  }
  
  template <class T>  
  bool compare(const T & exp_answer, const T & real_answer,
               const char * str = "result")
  {
    if (exp_answer != real_answer) {
      if (!m_end_of_line_printed) print_end_of_line();
      std::cout << "Expected " << str << ": " << exp_answer << std::endl
                << "Obtained " << str << ": " << real_answer << std::endl;
      m_end_of_line_printed = true;
      return false;
    } 
    return true;
  }

  template <class T>  
  bool compare_and_print(const T & exp_answer, const T & real_answer,
                         const char * str = "result")
  {
    bool result = compare(exp_answer, real_answer, str);
    print_result(result);
    return result;
  }
  
  template <class stream>
  bool read_point(stream & is, Point_2 & p);

  template <class stream>
  bool read_xcurve(stream & is, X_monotone_curve_2 & cv);
  
  template <class stream>
  bool read_curve(stream & is, Curve_2 & cv);
  
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
Traits_test<T_Traits>::Traits_test(int argc, char * argv[]) :
  m_end_of_line_printed(true)
{
  typedef T_Traits Traits;
  
  if (argc < 2 || argc >= 3) {
    std::cout << "Usage: " << argv[0] << " test_data_file" << std::endl;
    m_end_of_line_printed = true;
  } else
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
  m_filename.clear();
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
  std::ifstream is(m_filename.c_str());
  if (!is.is_open()) {
    std::cerr << "Error opening file " << m_filename.c_str() << std::endl;
    m_end_of_line_printed = true;
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
  unsigned int i;

  // Read points:
  unsigned int n_points;
  skip_comments(is, one_line);
  std::istringstream str_stream(one_line, std::istringstream::in);
  str_stream >> n_points;
  str_stream.clear();
  if (n_points > 0) {
    m_points.resize(n_points);
    for (i = 0; i < n_points; ++i) {
      skip_comments(is, one_line);
      str_stream.str(one_line);
      if (!read_point(str_stream, m_points[i])) {
        std::cerr << "Error reading point!" << std::endl;
        m_end_of_line_printed = true;
        return false;
      }
      str_stream.clear();
    }
  }

  // Read x-monotone curves
  unsigned int n_xcurves;
  skip_comments(is, one_line);
  str_stream.str(one_line);
  str_stream >> n_xcurves;
  str_stream.clear();
  if (n_xcurves > 0) {
    m_xcurves.resize(n_xcurves);
    for (i = 0; i < n_xcurves; ++i) {
      skip_comments(is, one_line);
      str_stream.str(one_line);
      if (!read_xcurve(str_stream, m_xcurves[i])) {
        std::cerr << "Error reading x-monotone curve!" << std::endl;
        m_end_of_line_printed = true;
        return false;
      }
      str_stream.clear();
    }
  }

  // Read curves
  unsigned int n_curves;
  skip_comments(is, one_line);
  str_stream.str(one_line);
  str_stream >> n_curves;
  str_stream.clear();
  if (n_curves > 0) {
    m_curves.resize(n_curves);
    for (i = 0; i < n_curves; ++i) {
      skip_comments(is, one_line);
      str_stream.str(one_line);
      if (!read_curve(str_stream, m_curves[i])) {
        std::cerr << "Error reading curve!" << std::endl;
        m_end_of_line_printed = true;
        return false;
      }
      str_stream.clear();
    }
  }

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
  m_end_of_line_printed = true;
  char one_line[128];
  char buff[128];
  while (!is.eof()) {
    skip_comments(is, one_line);
    std::istringstream str_stream(one_line, std::istringstream::in);
    buff[0] = '\0';
    str_stream.getline(buff, 128, ' ');
    std::string str_command(buff);
    Wrapper_iter wi = m_wrappers.find(str_command);
    str_stream.clear();
    if (wi == m_wrappers.end()) continue;
    Wrapper wrapper = (*wi).second;
    test_result &= (this->*wrapper)(str_stream);
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
  for (; *str != '\0'; ++str)
    if (*str != ' ') result += *str;
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

/*! Test Compare_x_2
 */
template <class T_Traits>
bool Traits_test<T_Traits>::compare_x_wrapper(std::istringstream & str_stream)
{
  unsigned int id1, id2;
  str_stream >> id1 >> id2;
  unsigned int exp_answer = get_expected_enum(str_stream);
  unsigned int real_answer =
    m_traits.compare_x_2_object()(m_points[id1], m_points[id2]);
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
  std::cout << "Test: compare_xy( " << m_points[id1] << ", "
            << m_points[id2] << " ) ? " << exp_answer << "  ";
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
  std::cout << "Test: min_vertex( " << m_xcurves[id1] << " ) ? "
            << exp_answer << " ";
  return compare_and_print(exp_answer, real_answer);
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
  std::cout << "Test: max_vertex( " << m_xcurves[id1] << " ) ? "
            << exp_answer << " ";
  return compare_and_print(exp_answer, real_answer);
}

template <class T_Traits>
bool
Traits_test<T_Traits>::is_vertical_wrapper(std::istringstream & str_stream)
{
  unsigned int id;
  str_stream >> id;
  bool exp_answer = get_expected_boolean(str_stream);
  bool real_answer = m_traits.is_vertical_2_object()(m_xcurves[id]);
  std::cout << "Test: is_vertical( " << m_xcurves[id]
            << " ) ? " << exp_answer << " ";
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
  std::cout << "Test: compare_y_at_x( " << m_points[id1] << ","
            << m_xcurves[id2]
            << " ) ? " << exp_answer << " ";
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
  unsigned int id1, id2, id3;
  str_stream >> id1 >> id2 >> id3;
  unsigned int exp_answer = get_expected_enum(str_stream);
  std::cout << "Test: compare_y_at_x_left( " << m_xcurves[id1] << ","
            << m_xcurves[id2] << ", " << m_points[id3]
            << " ) ? "; 
  unsigned int real_answer =
    m_traits.compare_y_at_x_left_2_object()(m_xcurves[id1], m_xcurves[id2],
                                            m_points[id3]);
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
    m_traits.compare_y_at_x_right_2_object()(m_xcurves[id1], m_xcurves[id2],
                                             m_points[id3]);
  std::cout << "Test: compare_y_at_x_right( " << m_xcurves[id1] << ","
            << m_xcurves[id2] << ", " << m_points[id3]
            << " ) ? " << exp_answer << " ";
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
  std::cout << "Test: equal( " << m_points[id1] << ", "
            << m_points[id2] << " ) ? " << exp_answer << " ";
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
  std::cout << "Test: equal( " << m_xcurves[id1] << ", "
            << m_xcurves[id2] << " ) ? " << exp_answer << " ";
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
  std::cout << "Test: make_x_monotone( " << m_curves[id]
            << " ) ? " << " ";

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
        m_end_of_line_printed = true;
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
      m_end_of_line_printed = true;
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
        m_end_of_line_printed = true;
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
      m_end_of_line_printed = true;
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
  std::cout << "Test: split( " << m_xcurves[id1] << "," << m_points[id2]
            << " ) ? " << " ";

  Equal_2 equal = m_traits.equal_2_object();
  str_stream >> id1 >> id2;                   // ... curve respectively

  if (!equal(m_xcurves[id1], cv1) || !equal(m_xcurves[id2], cv2)) {
    std::cerr << "Expected x-monotone curves: "
              << m_xcurves[id1] << ", " << m_xcurves[id2] << std::endl
              << "Obtained x-monotone curve: "
              << cv1 << "," << cv2 << std::endl;
    m_end_of_line_printed = true;
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
  unsigned int id1, id2;
  str_stream >> id1 >> id2;
  bool exp_answer = get_expected_boolean(str_stream);
  bool real_answer = m_traits.are_mergeable_2_object()(m_xcurves[id1],
                                                       m_xcurves[id2]);
  std::cout << "Test: are_mergeable( " << m_xcurves[id1] << ", "
            << m_xcurves[id2] << " ) ? " << exp_answer << " ";
  return compare_and_print(exp_answer, real_answer);
}

/*! Tests Merge_2.
 * Merge two given x-monotone curves into a single curve.
 */
template <class T_Traits>
bool Traits_test<T_Traits>::merge_wrapper(std::istringstream & str_stream)
{
  typedef T_Traits                              Traits;
  typedef typename Traits::X_monotone_curve_2   X_monotone_curve_2;
  typedef typename Traits::Equal_2              Equal_2;

  unsigned int id1, id2, id;
  str_stream >> id1 >> id2 >> id;
  X_monotone_curve_2 cv;
  m_traits.merge_2_object()(m_xcurves[id1], m_xcurves[id2], cv);
  std::cout << "Test: merge( " << m_xcurves[id1] << ", "
            << m_xcurves[id2] << " ) ? " << m_xcurves[id] << " ";
  Equal_2 equal = m_traits.equal_2_object();
  if (!equal(m_xcurves[id], cv)) {
    std::cerr << "Expected x-monotone curve: " << m_xcurves[id] << std::endl
              << "Obtained x-monotone curve: " << cv << std::endl;
    m_end_of_line_printed = true;
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
Traits_test<T_Traits>::approximate_wrapper(std::istringstream & str_stream)
{
  return false;
}

/*! tests Construct_x_monotone_curve_2.
 * Return an x-monotone curve connecting the two given endpoints.
 */
template <class T_Traits>
bool Traits_test<T_Traits>::
construct_x_monotone_curve_wrapper(std::istringstream & str_stream)
{
  return false;
}

template <class T_Traits>
template <class stream>
bool
Traits_test<T_Traits>::read_point(stream & is,
                                  typename T_Traits::Point_2 & p)
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
                                  typename T_Traits::X_monotone_curve_2 & cv)
{
  Basic_number_type x1, y1, x2, y2;
  is >> x1 >> y1 >> x2 >> y2;
  Point_2 p1(x1, y1);
  Point_2 p2(x2, y2);
  CGAL_assertion(p1 != p2);
  cv = typename T_Traits::X_monotone_curve_2(p1, p2);
  return true;
}

template <class T_Traits>
template <class stream>
bool
Traits_test<T_Traits>::read_curve(stream & is,
                                  typename T_Traits::Curve_2 & cv)
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
