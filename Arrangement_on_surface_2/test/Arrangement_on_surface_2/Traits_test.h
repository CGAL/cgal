#ifndef CGAL_TRAITS_TEST_H
#define CGAL_TRAITS_TEST_H

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
#include <CGAL/Arrangement_2/Arr_traits_adaptor_2_dispatching.h>

#include "Traits_base_test.h"

/*! Traits test */
template <typename T_Traits>
class Traits_test : public Traits_base_test<T_Traits> {
private:
  typedef T_Traits                                      Traits;
  typedef Traits_base_test<Traits>                      Base;
  typedef typename Base::Enum_type                      Enum_type;

  typedef typename Traits::Point_2                      Point_2;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename Traits::Curve_2                      Curve_2;
  
  /*! A map between (strings) commands and (member functions) operations */
  typedef bool(Traits_test::* Wrapper)(std::istringstream&);
  typedef std::map<std::string, Wrapper>                Wrapper_map;
  typedef typename Wrapper_map::iterator                Wrapper_iter;
  Wrapper_map m_wrappers;

  virtual bool exec(std::istringstream& str_stream,
                    const std::string& str_command,
                    bool& result)
  {
    Wrapper_iter wi = m_wrappers.find(str_command);
    str_stream.clear();
    if (wi == m_wrappers.end()) return true;
    Wrapper wrapper = (*wi).second;
    result = (this->*wrapper)(str_stream);
    return false;
  }

  //@{

  // traits

  /*! Test Compare_x_2
   */
  bool compare_x_wrapper(std::istringstream&);

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

  // /////////////////////////////////////////////////////////////////////////
  // boundary-specific functors:
  
  // -------------------------------------------------------------------------
  // left-right
  
  /*! Test Parameter_space_in_x_2
   */
  bool parameter_space_in_x_wrapper(std::istringstream&);
  bool parameter_space_in_x_wrapper_imp(std::istringstream&,
                                        CGAL::Arr_use_dummy_tag);
  bool parameter_space_in_x_wrapper_imp(std::istringstream&,
                                        CGAL::Arr_use_traits_tag);

  /*! Test Compare_y_near_boundary_2
   */
  bool compare_y_near_boundary_wrapper(std::istringstream&);
  bool compare_y_near_boundary_wrapper_imp(std::istringstream&,
                                           CGAL::Arr_use_dummy_tag);
  bool compare_y_near_boundary_wrapper_imp(std::istringstream&,
                                           CGAL::Arr_use_traits_tag);

  // TODO Is_on_y_identification_2
  
  
  // -------------------------------------------------------------------------
  // bottom-right
  
  /*! Test Parameter_space_in_y_2
   */
  bool parameter_space_in_y_wrapper(std::istringstream&);
  bool parameter_space_in_y_wrapper_imp(std::istringstream&,
                                        CGAL::Arr_use_dummy_tag);
  bool parameter_space_in_y_wrapper_imp(std::istringstream&,
                                        CGAL::Arr_use_traits_tag);

  /*! Test Compare_x_near_limit_2
   */
  bool compare_x_near_limit_wrapper(std::istringstream&);
  bool compare_x_near_limit_wrapper_imp(std::istringstream&,
                                        CGAL::Arr_use_dummy_tag);
  bool compare_x_near_limit_wrapper_imp(std::istringstream&,
                                        CGAL::Arr_use_traits_tag);

  /*! Test Compare_x_at_limit_2
   */
  bool compare_x_at_limit_wrapper(std::istringstream&);
  bool compare_x_at_limit_wrapper_imp(std::istringstream&,
                                      CGAL::Arr_use_dummy_tag);
  bool compare_x_at_limit_wrapper_imp(std::istringstream&,
                                      CGAL::Arr_use_traits_tag);

  /*! Test Compare_x_near_boundary_2
   */
  bool compare_x_near_boundary_wrapper(std::istringstream&);
  bool compare_x_near_boundary_wrapper_imp(std::istringstream&,
                                           CGAL::Arr_use_dummy_tag);
  bool compare_x_near_boundary_wrapper_imp(std::istringstream&,
                                           CGAL::Arr_use_traits_tag);

  /*! Test Compare_x_on_boundary_2
   */
  bool compare_x_on_boundary_wrapper(std::istringstream&);
  bool compare_x_on_boundary_wrapper_imp(std::istringstream&,
                                         CGAL::Arr_use_dummy_tag);
  bool compare_x_on_boundary_wrapper_imp(std::istringstream&,
                                         CGAL::Arr_use_traits_tag);
  
  // TODO Is_on_x_identification_2
  
  //@}

public:
  /*! Constructor */
  Traits_test(int argc, char * argv[]);

  /*! Destructor */
  ~Traits_test();
};

/*!
 * Constructor. 
 * Accepts test data file name.
 */
template <typename T_Traits>
Traits_test<T_Traits>::Traits_test(int argc, char * argv[]) :
  Traits_base_test<T_Traits>(argc, argv)
{
  typedef T_Traits Traits;
  
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

  // left-right
  
  m_wrappers[std::string("parameter_space_x")] =
    &Traits_test<Traits>::parameter_space_in_x_wrapper;
  m_wrappers[std::string("compare_y_near_boundary")] =
    &Traits_test<Traits>::compare_y_near_boundary_wrapper;

  // TODO Is_on_y_identification_2
  
  // bottom-top
  
  m_wrappers[std::string("parameter_space_y")] =
    &Traits_test<Traits>::parameter_space_in_y_wrapper;
  m_wrappers[std::string("compare_x_near_limit")] =
    &Traits_test<Traits>::compare_x_near_limit_wrapper;
  m_wrappers[std::string("compare_x_at_limit")] =
    &Traits_test<Traits>::compare_x_at_limit_wrapper;

  m_wrappers[std::string("compare_x_near_boundary")] =
    &Traits_test<Traits>::compare_x_near_boundary_wrapper;
  m_wrappers[std::string("compare_x_on_boundary")] =
    &Traits_test<Traits>::compare_x_on_boundary_wrapper;
  
  // TODO Is_on_x_identification_2
}

/*!
 * Destructor. 
 */
template <typename T_Traits>
Traits_test<T_Traits>::~Traits_test() {}


/*! Test Compare_x_2
 */
template <typename T_Traits>
bool Traits_test<T_Traits>::compare_x_wrapper(std::istringstream& str_stream)
{
  unsigned int id1, id2;
  str_stream >> id1 >> id2;
  unsigned int exp_answer = this->get_expected_enum(str_stream);
  std::cout << "Test: compare_x( " << this->m_points[id1] << ", "
            << this->m_points[id2] << " ) ? " << exp_answer << " ";

  unsigned int real_answer =
    this->m_traits.compare_x_2_object()(this->m_points[id1],
                                        this->m_points[id2]);
  return this->compare(exp_answer, real_answer);
}

/*! Test Compare_xy_2
 */
template <typename T_Traits>
bool Traits_test<T_Traits>::compare_xy_wrapper(std::istringstream& str_stream)
{
  unsigned int id1, id2;
  str_stream >> id1 >> id2;
  unsigned int exp_answer = this->get_expected_enum(str_stream);
  std::cout << "Test: compare_xy( " << this->m_points[id1] << ", "
            << this->m_points[id2] << " ) ? " << exp_answer << " ";

  unsigned int real_answer =
    this->m_traits.compare_xy_2_object()(this->m_points[id1],
                                         this->m_points[id2]);
  return this->compare(exp_answer, real_answer);
}

/*! Tests Construct_min_vertex_2.
 * Degenerate case: vertical curve.
 */
template <typename T_Traits>
bool Traits_test<T_Traits>::min_vertex_wrapper(std::istringstream& str_stream)
{
  unsigned int id1, id2;
  str_stream >> id1 >> id2;
  Point_2& exp_answer = this->m_points[id2];
  std::cout << "Test: min_vertex( " << this->m_xcurves[id1] << " ) ? "
            << exp_answer << " ";

  Point_2 real_answer =
    this->m_traits.construct_min_vertex_2_object()(this->m_xcurves[id1]);
  return this->compare_points(exp_answer, real_answer);
}

/*! Tests Construct_max_vertex_2.
 * Degenerate case: vertical curve.
 */
template <typename T_Traits>
bool Traits_test<T_Traits>::max_vertex_wrapper(std::istringstream& str_stream)
{
  unsigned int id1, id2;
  str_stream >> id1 >> id2;
  Point_2& exp_answer = this->m_points[id2];
  std::cout << "Test: max_vertex( " << this->m_xcurves[id1] << " ) ? "
            << exp_answer << " ";

  Point_2 real_answer =
    this->m_traits.construct_max_vertex_2_object()(this->m_xcurves[id1]);
  return this->compare_points(exp_answer, real_answer);
}

template <typename T_Traits>
bool
Traits_test<T_Traits>::is_vertical_wrapper(std::istringstream& str_stream)
{
  unsigned int id;
  str_stream >> id;
  bool exp_answer = this->get_expected_boolean(str_stream);
  std::cout << "Test: is_vertical( " << this->m_xcurves[id] << " ) ? "
            << exp_answer << " ";
  
  bool real_answer = this->m_traits.is_vertical_2_object()(this->m_xcurves[id]);
  return this->compare(exp_answer, real_answer);
}

/*! Tests Compare_y_at_x_2.
 * Return the location of the given point with respect to the input curve.
 * Degenerate cases: The point is an endpoint of the curve.
 *                   The curve is vertical.
 */
template <typename T_Traits>
bool
Traits_test<T_Traits>::compare_y_at_x_wrapper(std::istringstream& str_stream)
{
  unsigned int id1, id2;
  str_stream >> id1 >> id2;
  unsigned int exp_answer = this->get_expected_enum(str_stream);
  std::cout << "Test: compare_y_at_x( " << this->m_points[id1] << ","
            << this->m_xcurves[id2] << " ) ? " << exp_answer << " ";

  unsigned int real_answer =
    this->m_traits.compare_y_at_x_2_object()(this->m_points[id1],
                                             this->m_xcurves[id2]);
  return this->compare(exp_answer, real_answer);
}

/*! Tests Compare_y_at_x_left_2.
 * Compare the y value of two x-monotone curves immediately to the left
 * of their intersection point.
 * Degenerate cases: The curves coincide.
 *                   The curves coincide and vertical.
 *                   One of the curves is vertical.
 */
template <typename T_Traits>
bool
Traits_test<T_Traits>::
compare_y_at_x_left_wrapper(std::istringstream& str_stream)
{
  typedef typename T_Traits::Has_left_category          Has_left_category;
  return compare_y_at_x_left_wrapper_imp(str_stream, Has_left_category());
}

template <typename T_Traits>
bool
Traits_test<T_Traits>::
compare_y_at_x_left_wrapper_imp(std::istringstream&, CGAL::Tag_false)
{
  CGAL_error();
  return false;
}

template <typename T_Traits>
bool
Traits_test<T_Traits>::
compare_y_at_x_left_wrapper_imp(std::istringstream& str_stream, CGAL::Tag_true)
{
  unsigned int id1, id2, id3;
  str_stream >> id1 >> id2 >> id3;
  unsigned int exp_answer = this->get_expected_enum(str_stream);
  std::cout << "Test: compare_y_at_x_left( " << this->m_xcurves[id1] << ","
            << this->m_xcurves[id2] << ", " << this->m_points[id3] << " ) ? "
            << exp_answer << " ";

  unsigned int real_answer = this->m_traits.compare_y_at_x_left_2_object()
    (this->m_xcurves[id1], this->m_xcurves[id2], this->m_points[id3]);
  return this->compare(exp_answer, real_answer);
}

/*! Tests Compare_y_at_x_right_2.
 * Compare the y value of two x-monotone curves immediately to the right
 * of their intersection point.
 * Degenerate cases: The curves coincide.
 *                   The curves coincide and vertical.
 *                   One of the curves is vertical.
 */
template <typename T_Traits>
bool
Traits_test<T_Traits>::
compare_y_at_x_right_wrapper(std::istringstream& str_stream)
{
  unsigned int id1, id2, id3;
  str_stream >> id1 >> id2 >> id3;
  unsigned int exp_answer = this->get_expected_enum(str_stream);
  std::cout << "Test: compare_y_at_x_right( " << this->m_xcurves[id1] << ","
            << this->m_xcurves[id2] << ", " << this->m_points[id3] << " ) ? "
            << exp_answer << " ";

  unsigned int real_answer =
    this->m_traits.compare_y_at_x_right_2_object()
    (this->m_xcurves[id1], this->m_xcurves[id2], this->m_points[id3]);
  return this->compare(exp_answer, real_answer);
}

/*! Tests Equal_2::operator()(Point_2, Point_2).
 * Check whether two points are the same.
 */
template <typename T_Traits>
bool
Traits_test<T_Traits>::equal_points_wrapper(std::istringstream& str_stream)
{
  unsigned int id1, id2;
  str_stream >> id1 >> id2;
  bool exp_answer = this->get_expected_boolean(str_stream);
  std::cout << "Test: equal( " << this->m_points[id1] << ", "
            << this->m_points[id2] << " ) ? " << exp_answer << " ";

  bool real_answer = this->m_traits.equal_2_object()(this->m_points[id1], this->m_points[id2]);
  return this->compare(exp_answer, real_answer);
}

/*! Tests Equal_2::operator()(X_monotone_curve_2, X_monotone_curve_2).
 * Check whether two x-monotone curves are the same.
 */
template <typename T_Traits>
bool
Traits_test<T_Traits>::equal_curves_wrapper(std::istringstream& str_stream)
{
  unsigned int id1, id2;
  str_stream >> id1 >> id2;
  bool exp_answer = this->get_expected_boolean(str_stream);
  std::cout << "Test: equal( " << this->m_xcurves[id1] << ", "
            << this->m_xcurves[id2] << " ) ? " << exp_answer << " ";

  bool real_answer = this->m_traits.equal_2_object()(this->m_xcurves[id1], this->m_xcurves[id2]);
  return this->compare(exp_answer, real_answer);
}

/*! Tests Make_x_monotone_2.
 * Cut the given curve into x-monotone subcurves and insert them into the
 * given output iterator.
 * Degenerate cases for polylines: The first segment is vertical. The last
 * segment is vertical. Both firt and last are vertical. An internal segment
 * is vertical.
 */
template <typename T_Traits>
bool
Traits_test<T_Traits>::make_x_monotone_wrapper(std::istringstream& str_stream)
{
  typedef T_Traits                              Traits;
  typedef typename Traits::Point_2              Point_2;
  typedef typename Traits::X_monotone_curve_2   X_monotone_curve_2;
  typedef typename Traits::Curve_2              Curve_2;
  
  unsigned int id;
  str_stream >> id;
  std::cout << "Test: make_x_monotone( " << this->m_curves[id]
            << " ) ? ";
  std::vector<CGAL::Object> object_vec;
  this->m_traits.make_x_monotone_2_object()(this->m_curves[id],
                                      std::back_inserter(object_vec));
  unsigned int num;
  str_stream >> num;
  if (!this->compare(num, (unsigned int)object_vec.size(), "size"))
    return false;

  for (unsigned int i = 0; i < num; ++i) {
    unsigned int type;                  // 0 - point, 1 - x-monotone curve
    str_stream >> type;

    unsigned int id;                    // The id of the point or x-monotone
    str_stream >> id;                   // ... curve respectively

    unsigned int exp_type = 1;
    const X_monotone_curve_2 * xcv_ptr;
    xcv_ptr = CGAL::object_cast<X_monotone_curve_2> (&(object_vec[i]));
    if (xcv_ptr != NULL) {
      if (!this->compare(type, exp_type, "type")) return false;

      if (!this->compare_curves(this->m_xcurves[id], *xcv_ptr)) return false;
      continue;
    }

    exp_type = 0;
    const Point_2 * pt_ptr;
    pt_ptr = CGAL::object_cast<Point_2> (&(object_vec[i]));
    CGAL_assertion (pt_ptr != NULL);
    if (!this->compare(type, exp_type, "type")) return false;

    if (!this->compare_points(this->m_points[id], *pt_ptr)) return false;
  }
  object_vec.clear();
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
template <typename T_Traits>
bool Traits_test<T_Traits>::intersect_wrapper(std::istringstream& str_stream)
{
  typedef T_Traits                              Traits;
  typedef typename Traits::Point_2              Point_2;
  typedef typename Traits::X_monotone_curve_2   X_monotone_curve_2;
  typedef typename Traits::Curve_2              Curve_2;

  unsigned int id1, id2;
  str_stream >> id1 >> id2;
  std::vector<CGAL::Object> object_vec;
  this->m_traits.intersect_2_object()(this->m_xcurves[id1],
                                      this->m_xcurves[id2],
                                      std::back_inserter(object_vec));
  std::cout << "Test: intersect( " << this->m_xcurves[id1] << ","
            << this->m_xcurves[id2] << " ) ? ";
  std::size_t num;
  str_stream >> num;
  if (!this->compare(num, object_vec.size(), "size")) return false;

  for (unsigned int i = 0; i < num; ++i) {
    unsigned int type;                  // 0 - point, 1 - x-monotone curve
    str_stream >> type;
    unsigned int id;                    // The id of the point or x-monotone
    str_stream >> id;                   // ... curve respectively
    unsigned int multiplicity;
    if (type == 0) str_stream >> multiplicity;
    unsigned int exp_type = 1;
    const X_monotone_curve_2 * xcv_ptr =
      CGAL::object_cast<X_monotone_curve_2> (&(object_vec[i]));
    if (xcv_ptr != NULL) {
      if (!this->compare(type, exp_type, "type")) return false;

      if (!this->compare_curves(this->m_xcurves[id], *xcv_ptr)) return false;
      continue;
    }

    exp_type = 0;
    typedef std::pair<Point_2,unsigned int> Point_2_pair;
    const Point_2_pair * pt_pair_ptr =
      CGAL::object_cast<Point_2_pair> (&(object_vec[i]));
    CGAL_assertion(pt_pair_ptr != NULL);
    if (!this->compare(type, exp_type, "type")) return false;
    if (!this->compare_points(this->m_points[id], (*pt_pair_ptr).first))
      return false;
    if (!this->compare(multiplicity, (*pt_pair_ptr).second, "multiplicity"))
      return false;
  }
  object_vec.clear();
  return true;
}

/*! Tests Split_2.
 * Split a given x-monotone curve at a given point into two sub-curves.
 * Degenerate cases for polylines: the point and a polyline internal point
 * coincides.
 */
template <typename T_Traits>
bool Traits_test<T_Traits>::split_wrapper(std::istringstream& str_stream)
{
  typedef T_Traits                              Traits;
  typedef typename Traits::Point_2              Point_2;
  typedef typename Traits::X_monotone_curve_2   X_monotone_curve_2;
  typedef typename Traits::Curve_2              Curve_2;

  unsigned int id1, id2, id3, id4;
  str_stream >> id1 >> id2 >> id3 >> id4;
  X_monotone_curve_2 cv1, cv2;
  std::cout << "Test: split( " << this->m_xcurves[id1] << ","
            << this->m_points[id2] << " ) ? ";

  this->m_traits.split_2_object()(this->m_xcurves[id1], this->m_points[id2],
                                  cv1, cv2);
  return this->compare_curves(this->m_xcurves[id3], cv1) &&
    this->compare_curves(this->m_xcurves[id4], cv2);
}

/*! Tests Are_mergeable_2.
 * Check whether it is possible to merge two given x-monotone curves.
 */
template <typename T_Traits>
bool
Traits_test<T_Traits>::are_mergeable_wrapper(std::istringstream& str_stream)
{
  typedef typename T_Traits::Has_merge_category         Has_merge_category;
  return are_mergeable_wrapper_imp(str_stream, Has_merge_category());
}

template <typename T_Traits>
bool
Traits_test<T_Traits>::
are_mergeable_wrapper_imp(std::istringstream&, CGAL::Tag_false)
{
  CGAL_error();
  return false;
}

template <typename T_Traits>
bool Traits_test<T_Traits>::
are_mergeable_wrapper_imp(std::istringstream& str_stream, CGAL::Tag_true)
{
  unsigned int id1, id2;
  str_stream >> id1 >> id2;
  bool exp_answer = this->get_expected_boolean(str_stream);
  std::cout << "Test: are_mergeable( " << this->m_xcurves[id1] << ", "
            << this->m_xcurves[id2] << " ) ? " << exp_answer << " ";

  bool real_answer =
    this->m_traits.are_mergeable_2_object()(this->m_xcurves[id1],
                                            this->m_xcurves[id2]);
  return this->compare(exp_answer, real_answer);
}

/*! Tests Merge_2.
 * Merge two given x-monotone curves into a single curve.
 */
template <typename T_Traits>
bool Traits_test<T_Traits>::merge_wrapper(std::istringstream& str_stream)
{
  typedef typename T_Traits::Has_merge_category         Has_merge_category;
  return merge_wrapper_imp(str_stream, Has_merge_category());
}

template <typename T_Traits>
bool Traits_test<T_Traits>::merge_wrapper_imp(std::istringstream&,
                                              CGAL::Tag_false)
{
  CGAL_error();
  return false;
}

template <typename T_Traits>
bool Traits_test<T_Traits>::merge_wrapper_imp(std::istringstream& str_stream,
                                              CGAL::Tag_true)
{
  typedef T_Traits                              Traits;
  typedef typename Traits::X_monotone_curve_2   X_monotone_curve_2;

  unsigned int id1, id2, id;
  str_stream >> id1 >> id2 >> id;
  X_monotone_curve_2 cv;
  std::cout << "Test: merge( " << this->m_xcurves[id1] << ", "
            << this->m_xcurves[id2] << " ) ? " << this->m_xcurves[id] << " ";

  this->m_traits.merge_2_object()(this->m_xcurves[id1], this->m_xcurves[id2], cv);
  return this->compare_curves(this->m_xcurves[id], cv);
}

/*! Tests Approximate_2.
 * Return an approximation of a point coordinate.
 */
template <typename T_Traits>
bool
Traits_test<T_Traits>::approximate_wrapper(std::istringstream& )
{
  return false;
}

/*! tests Construct_x_monotone_curve_2.
 * Return an x-monotone curve connecting the two given endpoints.
 */
template <typename T_Traits>
bool Traits_test<T_Traits>::
construct_x_monotone_curve_wrapper(std::istringstream& )
{
  return false;
}

// ///////////////////////////////////////////////////////////////////////////
// boundary-specific functors


// ---------------------------------------------------------------------------
// left-right


/*! Test Parameter_space_in_x_2
*/
template <typename T_Traits>
bool Traits_test<T_Traits>::
parameter_space_in_x_wrapper(std::istringstream& str_stream)
{
  typedef typename CGAL::internal::Arr_complete_left_side_tag< T_Traits >::Tag
    Left_side_category;
  typedef typename CGAL::internal::Arr_complete_right_side_tag< T_Traits >::Tag
    Right_side_category;
  typedef CGAL::internal::Arr_left_right_implementation_dispatch
    <Left_side_category, Right_side_category>           LR;
  typedef typename LR::Parameter_space_in_x_2_curve_end_tag     Psx_tag;
  return parameter_space_in_x_wrapper_imp(str_stream, Psx_tag());
}

template <typename T_Traits>
bool Traits_test<T_Traits>::
parameter_space_in_x_wrapper_imp(std::istringstream&, CGAL::Arr_use_dummy_tag)
{
  CGAL_error();
  return false;
}

template <typename T_Traits>
bool Traits_test<T_Traits>::
parameter_space_in_x_wrapper_imp(std::istringstream& str_stream,
                                 CGAL::Arr_use_traits_tag)
{
  CGAL::Arr_parameter_space exp_answer, real_answer;
  unsigned int id;
  str_stream >> id;
  std::pair<Enum_type, unsigned int> next_input =
    this->get_next_input(str_stream);
  bool curves_op = (next_input.first == Base::CURVE_END);
  CGAL::Arr_curve_end cv_end = CGAL::ARR_MIN_END;
  if (curves_op) {
    cv_end = static_cast<CGAL::Arr_curve_end>(next_input.second);
    std::cout << "Test: parameter_space_x( " << this->m_xcurves[id] << " , "
              << (cv_end == CGAL::ARR_MIN_END ? "MIN_END" : "MAX_END")
              << " ) ? ";
    next_input = this->get_next_input(str_stream);
    CGAL_assertion(next_input.first == Base::PARAMETER_SPACE);
  } else if (next_input.first == Base::PARAMETER_SPACE) {
    std::cout << "Test: parameter_space_in_x( " << this->m_points[id] << " ) ? ";
  } else {
    CGAL_error();
    return false;
  }
  exp_answer = static_cast<CGAL::Arr_parameter_space>(next_input.second);
  std::cout << (exp_answer == CGAL::ARR_LEFT_BOUNDARY ? "LEFT_BOUNDARY" :
               (exp_answer == CGAL::ARR_RIGHT_BOUNDARY ? "RIGHT_BOUNDARY" :
               (exp_answer == CGAL::ARR_BOTTOM_BOUNDARY ? "BOTTOM_BOUNDARY" :
               (exp_answer == CGAL::ARR_TOP_BOUNDARY ?
                "TOP_BOUNDARY":"INTERIOR")))) << " ";

  real_answer = (curves_op) ?
    this->m_traits.parameter_space_in_x_2_object()(this->m_xcurves[id],
                                                   cv_end) :
    this->m_traits.parameter_space_in_x_2_object()(this->m_points[id]);
  return this->compare(exp_answer, real_answer);
}

/*! Test Compare_y_near_boundary_2
 */
template <typename T_Traits>
bool Traits_test<T_Traits>::
compare_y_near_boundary_wrapper(std::istringstream& str_stream)
{
  typedef typename CGAL::internal::Arr_complete_left_side_tag< T_Traits >::Tag
    Left_side_category;
  typedef typename CGAL::internal::Arr_complete_right_side_tag< T_Traits >::Tag
    Right_side_category;
  typedef CGAL::internal::Arr_left_right_implementation_dispatch
    <Left_side_category, Right_side_category>           LR;
  typedef typename LR::Compare_y_near_boundary_2_curve_ends_tag Cmp_tag;
  return compare_y_near_boundary_wrapper_imp(str_stream, Cmp_tag());
}

template <typename T_Traits>
bool Traits_test<T_Traits>::
compare_y_near_boundary_wrapper_imp(std::istringstream&, CGAL::Arr_use_dummy_tag)
{
  CGAL_error();
  return false;
}

template <typename T_Traits>
bool Traits_test<T_Traits>::
compare_y_near_boundary_wrapper_imp(std::istringstream& str_stream,
                                    CGAL::Arr_use_traits_tag)
{
  unsigned int id1, id2;
  str_stream >> id1 >> id2;
  std::pair<Enum_type, unsigned int> next_input =
    this->get_next_input(str_stream);
  CGAL_assertion(next_input.first == Base::CURVE_END);
  CGAL::Arr_curve_end cv_end =
    static_cast<CGAL::Arr_curve_end>(next_input.second);
  std::cout << "Test: compare_y_near_boundary( " << this->m_xcurves[id1] << " , "
            << this->m_xcurves[id2]<< " , "
            << (cv_end == CGAL::ARR_MIN_END ? "MIN_END" : "MAX_END") << " ) ? ";
  next_input = this->get_next_input(str_stream);
  CGAL_assertion(next_input.first == Base::SIGN);
  CGAL::Comparison_result exp_answer = 
    static_cast<CGAL::Comparison_result>(next_input.second);
  std::cout << (exp_answer == CGAL::SMALLER ? "SMALLER":
               (exp_answer == CGAL::LARGER ? "LARGER" : "EQUAL")) << " ";

  CGAL::Comparison_result real_answer = 
    this->m_traits.compare_y_near_boundary_2_object()
    (this->m_xcurves[id1], this->m_xcurves[id2], cv_end);
  return this->compare(exp_answer, real_answer);
}

// TODO Is_on_y_identification_2

// ---------------------------------------------------------------------------
// bottom-top

/*! Test Parameter_space_in_y_2
 */
template <typename T_Traits>
bool Traits_test<T_Traits>::
parameter_space_in_y_wrapper(std::istringstream& str_stream)
{
  typedef typename CGAL::internal::Arr_complete_bottom_side_tag<T_Traits>::Tag
    Bottom_side_category;
  typedef typename CGAL::internal::Arr_complete_top_side_tag<T_Traits>::Tag
    Top_side_category;
  typedef CGAL::internal::Arr_bottom_top_implementation_dispatch
    <Bottom_side_category, Top_side_category>           BT;
  typedef typename BT::Parameter_space_in_y_2_curve_end_tag     Psy_tag;
  return parameter_space_in_y_wrapper_imp(str_stream, Psy_tag());
}

template <typename T_Traits>
bool Traits_test<T_Traits>::
parameter_space_in_y_wrapper_imp(std::istringstream&, CGAL::Arr_use_dummy_tag)
{
  CGAL_error();
  return false;
}

template <typename T_Traits>
bool Traits_test<T_Traits>::
parameter_space_in_y_wrapper_imp(std::istringstream& str_stream,
                                 CGAL::Arr_use_traits_tag)
{
  unsigned int id;
  str_stream >> id;
  std::pair<Enum_type, unsigned int> next_input =
    this->get_next_input(str_stream);
  bool curves_op = (next_input.first == Base::CURVE_END);
  CGAL::Arr_curve_end cv_end = CGAL::ARR_MIN_END;
  if (curves_op) {
    cv_end = static_cast<CGAL::Arr_curve_end>(next_input.second);
    std::cout << "Test: parameter_space_y( " << this->m_xcurves[id] << " , "
              << (cv_end == CGAL::ARR_MIN_END ? "MIN_END" : "MAX_END")
              << " ) ? ";
    next_input = this->get_next_input(str_stream);
    CGAL_assertion(next_input.first == Base::PARAMETER_SPACE);
  } else if (next_input.first == Base::PARAMETER_SPACE) {
    std::cout << "Test: parameter_space_in_y( " << this->m_points[id] << " ) ? ";
  } else {
    CGAL_error();
    return false;
  }
  CGAL::Arr_parameter_space exp_answer =
    static_cast<CGAL::Arr_parameter_space>(next_input.second);
  std::cout << (exp_answer == CGAL::ARR_LEFT_BOUNDARY ? "LEFT_BOUNDARY":
               (exp_answer == CGAL::ARR_RIGHT_BOUNDARY ? "RIGHT_BOUNDARY":
               (exp_answer == CGAL::ARR_BOTTOM_BOUNDARY ? "BOTTOM_BOUNDARY":
               (exp_answer == CGAL::ARR_TOP_BOUNDARY ?
                "TOP_BOUNDARY" : "INTERIOR")))) << " ";

  CGAL::Arr_parameter_space real_answer = (curves_op) ?
    this->m_traits.parameter_space_in_y_2_object()(this->m_xcurves[id],
                                                   cv_end) :
    this->m_traits.parameter_space_in_y_2_object()(this->m_points[id]);
  return this->compare(exp_answer, real_answer);
}

/* Compare_x_near_limit_2
 */
template <typename T_Traits>
bool Traits_test<T_Traits>::
compare_x_near_limit_wrapper(std::istringstream& str_stream)
{
  typedef typename CGAL::internal::Arr_complete_bottom_side_tag<T_Traits>::Tag
    Bottom_side_category;
  typedef typename CGAL::internal::Arr_complete_top_side_tag<T_Traits>::Tag
    Top_side_category;
  typedef CGAL::internal::Arr_bottom_top_implementation_dispatch
    <Bottom_side_category, Top_side_category>           BT;
  typedef typename BT::Compare_x_near_limit_2_curve_ends_tag Cmp_tag;
  return compare_x_near_limit_wrapper_imp(str_stream, Cmp_tag());
}

template <typename T_Traits>
bool Traits_test<T_Traits>::
compare_x_near_limit_wrapper_imp(std::istringstream&,
                                 CGAL::Arr_use_dummy_tag)
{
  CGAL_error();
  return false;
}

template <typename T_Traits>
bool Traits_test<T_Traits>::
compare_x_near_limit_wrapper_imp(std::istringstream& str_stream,
                                 CGAL::Arr_use_traits_tag)
{
  std::cout << "Test: compare_x_near_limit( ";

  unsigned int id1, id2;
  // xcurve xcurve curve-end
  str_stream >> id1 >> id2;
  std::pair<Enum_type, unsigned int> next_input =
    this->get_next_input(str_stream);
  CGAL_assertion(next_input.first == Base::CURVE_END);
  CGAL::Arr_curve_end cv_end =
    static_cast<CGAL::Arr_curve_end>(next_input.second);
  
  std::cout << this->m_xcurves[id1] << ", "
            << this->m_xcurves[id2] << " , "
            << this->curve_end_str(cv_end) << " ) ? ";

  next_input = this->get_next_input(str_stream);
  CGAL_assertion(next_input.first == Base::SIGN);
  CGAL::Comparison_result exp_answer = 
    static_cast<CGAL::Comparison_result>(next_input.second);
  std::cout << (exp_answer == CGAL::SMALLER ? "SMALLER":
               (exp_answer == CGAL::LARGER ? "LARGER":"EQUAL")) << " ";
  
  CGAL::Comparison_result real_answer = 
    this->m_traits.compare_x_near_limit_2_object()
    (this->m_xcurves[id1], this->m_xcurves[id2], cv_end);
  return this->compare(exp_answer, real_answer);
}

/* Compare_x_at_limit
 */
template <typename T_Traits>
bool Traits_test<T_Traits>::
compare_x_at_limit_wrapper(std::istringstream& str_stream)
{
  typedef typename CGAL::internal::Arr_complete_bottom_side_tag<T_Traits>::Tag
    Bottom_side_category;
  typedef typename CGAL::internal::Arr_complete_top_side_tag<T_Traits>::Tag
    Top_side_category;
  typedef CGAL::internal::Arr_bottom_top_implementation_dispatch
    <Bottom_side_category, Top_side_category> BT;
  typedef typename BT::Compare_x_at_limit_2_curve_ends_tag             Cmp_tag1;
  typedef typename BT::Compare_x_at_limit_2_point_curve_end_tag        Cmp_tag2;
  typedef typename CGAL::internal::Or_traits<Cmp_tag1, Cmp_tag2>::type Cmp_tag;
  return compare_x_at_limit_wrapper_imp(str_stream, Cmp_tag());
}

template <typename T_Traits>
bool Traits_test<T_Traits>::
compare_x_at_limit_wrapper_imp(std::istringstream&,
                               CGAL::Arr_use_dummy_tag)
{
  CGAL_error();
  return false;
}

template <typename T_Traits>
bool Traits_test<T_Traits>::
compare_x_at_limit_wrapper_imp(std::istringstream& str_stream,
                               CGAL::Arr_use_traits_tag)
{
  std::cout << "Test: compare_x_at_limit( ";
  
  CGAL::Comparison_result real_answer;
  unsigned int id1, id2;
  // first argument must be a number (either a point or a xcurve)
  str_stream >> id1;
  std::pair<Enum_type, unsigned int> next_input =
    this->get_next_input(str_stream);
  // second argument can be number or text (either xcurve or curve_end)
  bool curves_op = next_input.first == Base::NUMBER;
  CGAL::Arr_curve_end cv_end1 = CGAL::ARR_MIN_END, cv_end2 = CGAL::ARR_MIN_END;

  if (curves_op) {
    id2 = static_cast<unsigned int>(next_input.second);
    next_input = this->get_next_input(str_stream);
    CGAL_assertion(next_input.first == Base::CURVE_END);
    cv_end1 = static_cast<CGAL::Arr_curve_end>(next_input.second);
      
    std::cout << this->m_points[id1] << " , "
              << this->m_xcurves[id2]<< " , "
              << this->curve_end_str(cv_end1) << " ) ? ";
  }
  else if (next_input.first == Base::CURVE_END) {
    cv_end1 = static_cast<CGAL::Arr_curve_end>(next_input.second);
    next_input = this->get_next_input(str_stream);
    CGAL_assertion(next_input.first == Base::NUMBER);
    id2 = static_cast<unsigned int>(next_input.second);
    next_input = this->get_next_input(str_stream);
    CGAL_assertion(next_input.first == Base::CURVE_END);
    cv_end2 = static_cast<CGAL::Arr_curve_end>(next_input.second);

    std::cout << this->m_xcurves[id1] << ", "
              << this->curve_end_str(cv_end1) << ", "
              << this->m_xcurves[id2] << " , "
              << this->curve_end_str(cv_end2) << " ) ? ";
  }
  else CGAL_error();

  next_input = this->get_next_input(str_stream);
  CGAL_assertion(next_input.first == Base::SIGN);
  CGAL::Comparison_result exp_answer = 
    static_cast<CGAL::Comparison_result>(next_input.second);
  std::cout << (exp_answer == CGAL::SMALLER ? "SMALLER":
               (exp_answer == CGAL::LARGER ? "LARGER":"EQUAL")) << " ";

  real_answer = (curves_op) ?
    this->m_traits.compare_x_at_limit_2_object()
    (this->m_points[id1], this->m_xcurves[id2], cv_end1) :
    this->m_traits.compare_x_at_limit_2_object()
    (this->m_xcurves[id1], cv_end1,this->m_xcurves[id2], cv_end2);
  return this->compare(exp_answer, real_answer);
}

/* Compare_x_near_boundary_2
 */
template <typename T_Traits>
bool Traits_test<T_Traits>::
compare_x_near_boundary_wrapper(std::istringstream& str_stream)
{
  typedef typename CGAL::internal::Arr_complete_bottom_side_tag<T_Traits>::Tag
    Bottom_side_category;
  typedef typename CGAL::internal::Arr_complete_top_side_tag<T_Traits>::Tag
    Top_side_category;
  typedef CGAL::internal::Arr_bottom_top_implementation_dispatch
    <Bottom_side_category, Top_side_category>           BT;
  typedef typename BT::Compare_x_near_boundary_2_curve_ends_tag Cmp_tag;
  return compare_x_near_boundary_wrapper_imp(str_stream, Cmp_tag());
}

template <typename T_Traits>
bool Traits_test<T_Traits>::
compare_x_near_boundary_wrapper_imp(std::istringstream&,
                                    CGAL::Arr_use_dummy_tag)
{
  CGAL_error();
  return false;
}

template <typename T_Traits>
bool Traits_test<T_Traits>::
compare_x_near_boundary_wrapper_imp(std::istringstream& str_stream,
                                    CGAL::Arr_use_traits_tag)
{
  std::cout << "Test: compare_x_near_boundary( ";

  unsigned int id1, id2;
  // xcurve xcurve curve-end
  str_stream >> id1 >> id2;
  std::pair<Enum_type, unsigned int> next_input =
    this->get_next_input(str_stream);
  CGAL_assertion(next_input.first == Base::CURVE_END);
  CGAL::Arr_curve_end cv_end =
    static_cast<CGAL::Arr_curve_end>(next_input.second);
  
  std::cout << this->m_xcurves[id1] << ", "
            << this->m_xcurves[id2] << " , "
            << this->curve_end_str(cv_end) << " ) ? ";

  next_input = this->get_next_input(str_stream);
  CGAL_assertion(next_input.first == Base::SIGN);
  CGAL::Comparison_result exp_answer = 
    static_cast<CGAL::Comparison_result>(next_input.second);
  std::cout << (exp_answer == CGAL::SMALLER ? "SMALLER":
               (exp_answer == CGAL::LARGER ? "LARGER":"EQUAL")) << " ";
  
  CGAL::Comparison_result real_answer = 
    this->m_traits.compare_x_near_boundary_2_object()
    (this->m_xcurves[id1], this->m_xcurves[id2], cv_end);    
  return this->compare(exp_answer, real_answer);
}

/* Compare_x_on_boundary
 */
template <typename T_Traits>
bool Traits_test<T_Traits>::
compare_x_on_boundary_wrapper(std::istringstream& str_stream)
{
  typedef typename CGAL::internal::Arr_complete_bottom_side_tag<T_Traits>::Tag
    Bottom_side_category;
  typedef typename CGAL::internal::Arr_complete_top_side_tag<T_Traits>::Tag
    Top_side_category;
  typedef CGAL::internal::Arr_bottom_top_implementation_dispatch
    <Bottom_side_category, Top_side_category> BT;
  typedef typename BT::Compare_x_on_boundary_2_points_tag              Cmp_tag1;
  typedef typename BT::Compare_x_on_boundary_2_point_curve_end_tag     Cmp_tag2;
  typedef typename BT::Compare_x_on_boundary_2_curve_ends_tag          Cmp_tag3;
  typedef typename CGAL::internal::Or_traits<Cmp_tag1, Cmp_tag2>::type Cmp_tag12;
  typedef typename CGAL::internal::Or_traits<Cmp_tag12, Cmp_tag3>::type Cmp_tag;
  return compare_x_on_boundary_wrapper_imp(str_stream, Cmp_tag());
}

template <typename T_Traits>
bool Traits_test<T_Traits>::
compare_x_on_boundary_wrapper_imp(std::istringstream&,
                                  CGAL::Arr_use_dummy_tag)
{
  CGAL_error();
  return false;
}

template <typename T_Traits>
bool Traits_test<T_Traits>::
compare_x_on_boundary_wrapper_imp(std::istringstream& str_stream,
                                  CGAL::Arr_use_traits_tag)
{
  // TBD: points
  std::cout << "Test: compare_x_on_boundary( ";
  
  CGAL::Comparison_result real_answer;
  unsigned int id1, id2;
  // first argument must be a number (either a point or a xcurve)
  str_stream >> id1;
  std::pair<Enum_type, unsigned int> next_input =
    this->get_next_input(str_stream);
  // second argument can be number or text (either xcurve or curve_end)
  bool curves_op = next_input.first == Base::NUMBER;
  CGAL::Arr_curve_end cv_end1 = CGAL::ARR_MIN_END, cv_end2 = CGAL::ARR_MIN_END;

  if (curves_op) {
    id2 = static_cast<unsigned int>(next_input.second);
    next_input = this->get_next_input(str_stream);
    CGAL_assertion(next_input.first == Base::CURVE_END);
    cv_end1 = static_cast<CGAL::Arr_curve_end>(next_input.second);
      
    std::cout << this->m_points[id1] << " , "
              << this->m_xcurves[id2]<< " , "
              << this->curve_end_str(cv_end1) << " ) ? ";
  }
  else if (next_input.first == Base::CURVE_END) {
    cv_end1 = static_cast<CGAL::Arr_curve_end>(next_input.second);
    next_input = this->get_next_input(str_stream);
    CGAL_assertion(next_input.first == Base::NUMBER);
    id2 = static_cast<unsigned int>(next_input.second);
    next_input = this->get_next_input(str_stream);
    CGAL_assertion(next_input.first == Base::CURVE_END);
    cv_end2 = static_cast<CGAL::Arr_curve_end>(next_input.second);

    std::cout << this->m_xcurves[id1] << ", "
              << this->curve_end_str(cv_end1) << ", "
              << this->m_xcurves[id2] << " , "
              << this->curve_end_str(cv_end2) << " ) ? ";
  }
  else CGAL_error();

  next_input = this->get_next_input(str_stream);
  CGAL_assertion(next_input.first == Base::SIGN);
  CGAL::Comparison_result exp_answer = 
    static_cast<CGAL::Comparison_result>(next_input.second);
  std::cout << (exp_answer == CGAL::SMALLER ? "SMALLER":
               (exp_answer == CGAL::LARGER ? "LARGER":"EQUAL")) << " ";

  real_answer = (curves_op) ?
    this->m_traits.compare_x_on_boundary_2_object()
    (this->m_points[id1], this->m_xcurves[id2], cv_end1) :
    this->m_traits.compare_x_on_boundary_2_object()
    (this->m_xcurves[id1], cv_end1,this->m_xcurves[id2], cv_end2);
  return this->compare(exp_answer, real_answer);
}

// TODO Is_on_x_identification_2

#endif
