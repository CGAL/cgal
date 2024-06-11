#ifndef CGAL_TRAITS_TEST_H
#define CGAL_TRAITS_TEST_H


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
// #include <cstdlib>

#include <boost/lexical_cast.hpp>

#include <CGAL/exceptions.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arrangement_2/Arr_traits_adaptor_2_dispatching.h>
#include <CGAL/use.h>

#include "Traits_base_test.h"

/*! Traits test */
template <typename GeomTraits>
class Traits_test : public Traits_base_test<GeomTraits> {
private:
  Traits_test<GeomTraits>& operator=(const Traits_test<GeomTraits>&);
  typedef GeomTraits                                    Traits;
  typedef Traits_base_test<Traits>                      Base;
  typedef typename Base::Enum_type                      Enum_type;

  typedef typename Traits::Point_2                      Point_2;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename Traits::Curve_2                      Curve_2;

  // some polycurve functors needs Segment and x-monotone segment to be defined
  // which are normally not found in other geom_traits.
#if TEST_GEOM_TRAITS == POLYCURVE_CONIC_GEOM_TRAITS ||          \
  TEST_GEOM_TRAITS == POLYCURVE_CIRCULAR_ARC_GEOM_TRAITS ||     \
  TEST_GEOM_TRAITS == POLYCURVE_BEZIER_GEOM_TRAITS

  typedef typename Traits::Subcurve_2                   Subcurve_2;
  typedef typename Traits::X_monotone_subcurve_2        X_monotone_subcurve_2;

#endif

  /*! A map between (strings) commands and (member functions) operations */
  typedef bool(Traits_test::* Wrapper)(std::istringstream&);
  typedef std::map<std::string, Wrapper>                Wrapper_map;
  typedef typename Wrapper_map::iterator                Wrapper_iter;
  Wrapper_map m_wrappers;

  virtual bool exec(std::istringstream& str_stream,
                    const std::string& str_command,
                    bool& result) {
    // str_stream is the input file object.
    // Get the appropriate functor. "str_command" consist of the appropriate
    // functor string.
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
   * segment is vertical. Both first and last are vertical. An internal segment
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

  /*! Fall threw implementation to handle the case where the
   * Approximate_2::operator()(const X_monotone_curve_2& xcv,
   *                           Approximate_number_type error,
   *                           OutputIterator oi,
   *                           bool l2r)
   * is not supported.
   */
  template <typename GT>
  bool approximate_wrapper_impl2(std::istringstream& is, long);

  /*! Specialized instance to handle the case where the
   * Approximate_2::operator()(const X_monotone_curve_2& xcv,
   *                           Approximate_number_type error,
   *                           OutputIterator oi,
   *                           bool l2r)
   * is supported.
   */
  template <typename GT,
            typename T =
            decltype(std::declval<GT>().approximate_2_object().operator()
                     (std::declval<typename GT::X_monotone_curve_2>(),
                      0,
                      (std::declval<int*>())))>
  bool approximate_wrapper_impl2(std::istringstream& is, int);

  /*! Fall threw implementation to handle the case where the type
   * Traits::Approximate_2 does not exist.
   */
  template <typename GT>
  bool approximate_wrapper_impl1(std::istringstream& is, long);

  /*! Specialized instance to handle the case where the type
   * Traits::Approximate_2 exists.
   */
  template <typename GT,
            typename T = decltype(std::declval<GT>().approximate_2_object())>
  bool approximate_wrapper_impl1(std::istringstream& is, int);

  /*!
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

  /*
  * Test Push_back
  */
  // some polycurve functors needs Segment and x-monotone segment to be defined
  // which are normally not found in other geom_traits.
  #if TEST_GEOM_TRAITS == POLYCURVE_CONIC_GEOM_TRAITS || \
      TEST_GEOM_TRAITS == POLYCURVE_CIRCULAR_ARC_GEOM_TRAITS || \
      TEST_GEOM_TRAITS == POLYCURVE_BEZIER_GEOM_TRAITS || \
      TEST_GEOM_TRAITS == POLYLINE_GEOM_TRAITS
  bool push_back_wrapper(std::istringstream& str_stream);
  bool push_front_wrapper(std::istringstream& str_stream);
  bool number_of_points_wrapper(std::istringstream& str_stream);
  bool compare_endpoints_xy_wrapper(std::istringstream& str_stream);
  bool construct_opposite_wrapper(std::istringstream& str_stream);
  bool trim_wrapper(std::istringstream& str_stream);
  #endif
  // TODO Is_on_x_identification_2

  //@}

public:
  /*! Constructor */
  Traits_test(const Traits& traits);

  /*! Destructor */
  ~Traits_test();
};

/*! Constructor.
 * Accepts test data file name.
 */
template <typename GeomTraits>
Traits_test<GeomTraits>::Traits_test(const GeomTraits& traits) : Base(traits) {
  using Traits = GeomTraits;

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

  m_wrappers[std::string("compare_x_near_boundary")] =
    &Traits_test<Traits>::compare_x_near_boundary_wrapper;
  m_wrappers[std::string("compare_x_on_boundary")] =
    &Traits_test<Traits>::compare_x_on_boundary_wrapper;

  // some polycurve functors needs Segment and x-monotone segment to be defined
  // which are normally not found in other geom_traits.
#if TEST_GEOM_TRAITS == POLYCURVE_CONIC_GEOM_TRAITS || \
      TEST_GEOM_TRAITS == POLYCURVE_CIRCULAR_ARC_GEOM_TRAITS || \
      TEST_GEOM_TRAITS == POLYCURVE_BEZIER_GEOM_TRAITS || \
      TEST_GEOM_TRAITS == POLYLINE_GEOM_TRAITS
  m_wrappers[std::string("push_back")] =
    &Traits_test<Traits>::push_back_wrapper;
  m_wrappers[std::string("push_front")] =
    &Traits_test<Traits>::push_front_wrapper;
  m_wrappers[std::string("number_of_points")] =
    &Traits_test<Traits>::number_of_points_wrapper;
  m_wrappers[std::string("compare_endpoints_xy")] =
    &Traits_test<Traits>::compare_endpoints_xy_wrapper;
  m_wrappers[std::string("construct_opposite")] =
    &Traits_test<Traits>::construct_opposite_wrapper;
  m_wrappers[std::string("trim")] =
    &Traits_test<Traits>::trim_wrapper;
#endif
  // TODO Is_on_x_identification_2
}

/*! Destructor.
 */
template <typename GeomTraits>
Traits_test<GeomTraits>::~Traits_test() {}

// some polycurve functors needs Segment and x-monotone segment to be defined
// which are normally not found in other geom_traits.
#if TEST_GEOM_TRAITS == POLYCURVE_CONIC_GEOM_TRAITS || \
    TEST_GEOM_TRAITS == POLYCURVE_CIRCULAR_ARC_GEOM_TRAITS || \
    TEST_GEOM_TRAITS == POLYCURVE_BEZIER_GEOM_TRAITS || \
    TEST_GEOM_TRAITS == POLYLINE_GEOM_TRAITS

template <typename GeomTraits>
bool Traits_test<GeomTraits>::trim_wrapper(std::istringstream& str_stream) {
  unsigned int x_curve_id, xcv_trimmed, src_id, tgt_id;

  // Read the ID's of the x-curve, source and target points
  // and the trimmed xcv.
  str_stream >> x_curve_id >> src_id >> tgt_id >> xcv_trimmed;

  //get the x-monotone curve
  X_monotone_curve_2 xcv = this->m_xcurves[x_curve_id];

  //get the trimmed curve for confirmation.
  X_monotone_curve_2 expected_xcv = this->m_xcurves[xcv_trimmed];

  //get the trimming source and target points
  Point_2 src = this->m_points[src_id];
  Point_2 tgt = this->m_points[tgt_id];

  std::cout << "Test: Trim ( " << xcv << " from "
            << src << " to " << tgt << " ) ?";

  X_monotone_curve_2 trimmed_xcv =
    this->m_geom_traits.trim_2_object()(xcv, src, tgt);

  if (!this->compare_curves(trimmed_xcv, expected_xcv)) { return false; }

  return true;
}

/* Test Push_back
 */
template <typename GeomTraits>
bool Traits_test<GeomTraits>::push_back_wrapper(std::istringstream& str_stream)
{
  //type: 0 for pushing a segment into curve.
  //      1 for pushing x-monotone segment into x-monotone curve.
  unsigned int type;
  str_stream >> type;

  // Ids of base curve/x-curve.
  unsigned int id1;
  str_stream >> id1;

  unsigned int segment_id;
  str_stream >> segment_id;

  //id of expected polycurve/x-monotone polycurve
  unsigned int expected_curve_id;
  str_stream >> expected_curve_id;

  if (type == 0) {
    /* THERE IS NO WAY AS OF NOW TO CHECK IF THE POLYCURVE (NON X-MONOTONE) IS
     * EQUAL. HENCE, UNTIL THAT COMPARISON IS NOT AVAILABLE IN THE
     * ARR_POLYCURVE_TRAITS, THIS TEST WILL PASS ONLY IF THE PRINTED RESULT
     * OF THE EXPECTED CURVE AND THE ACTUAL OBTAINED CURVE IS IDENTICAL.
     */
#if 0
    Curve_2 base_curve = this->m_curves[id1];
    Subcurve_2 segment = this->m_segments[segment_id];
    std::cout << "Test: push_back ( " << segment << " into "
              << base_curve << " ) ? ";
    this->m_geom_traits.push_back_2_object()( base_curve, segment );
    Curve_2 exp_curve = this->m_curves[expected_curve_id];
    std::stringstream sstr1, sstr2;
    sstr1 << std::cout << base_curve;
    sstr2 << std::cout << exp_curve;
    if (sstr1.str() != sstr2.str()) {
      std::cout << "Obtained result and expected result does not match"
                << std::endl;
      std::cout << std::endl << "Result obtained: " << sstr1.str() << std::endl;
      std::cout << std::endl << "Expected result: " << sstr2.str() << std::endl;
      return false;
    }
#endif
  }

  else if (type == 1) {
    X_monotone_curve_2 base_curve = this->m_xcurves[id1];
    X_monotone_subcurve_2 x_segment = this->m_xsegments[segment_id];

    std::cout << "Test: push_back ( "
              << x_segment << " into "
              << base_curve << " ) ? ";

    this->m_geom_traits.push_back_2_object()( base_curve, x_segment );

    X_monotone_curve_2 exp_curve = this->m_xcurves[expected_curve_id];

    if (!this->compare_curves(exp_curve, base_curve)) return false;
  }
  else {
    std::cout << "Incorrect type of operator. "
              << "Please refer to the descriptopn in the data file."
              << std::endl;
    return false;
  }

  return true;
}

/*
 * Test Push_front
 */
template <typename GeomTraits>
bool Traits_test<GeomTraits>::
push_front_wrapper(std::istringstream& str_stream) {
  //type: 0 for pushing a segment into curve.
  //      1 for pushing x-monotone segment into x-monotone curve.
  unsigned int type;
  str_stream >> type;

  // Ids of base curve/x-curve.
  unsigned int id1;
  str_stream >> id1;

  unsigned int segment_id;
  str_stream >> segment_id;

  //id of expected polycurve/x-monotone polycurve
  unsigned int expected_curve_id;
  str_stream >> expected_curve_id;

  if (type == 0) {
    /* THERE IS NO WAY AS OF NOW TO CHECK IF THE POLYCURVE (NON X-MONOTONE) IS
     * EQUAL. HENCE, UNTIL THAT COMPARISON IS NOT AVAILABLE IN THE
     * ARR_POLYCURVE_TRAITS, THIS TEST WILL PASS ONLY IF THE PRINTED RESULT
     * OF THE EXPECTED CURVE AND THE ACTUAL OBTAINED CURVE IS IDENTICAL.
     */
#if 0
    Curve_2 base_curve = this->m_curves[id1];
    Subcurve_2 segment = this->m_segments[segment_id];
    std::cout << "Test: push_front ( " << segment << "into "
              << base_curve << " ) ? ";
    this->m_geom_traits.push_front_2_object()( base_curve, segment );
    Curve_2 exp_curve = this->m_curves[expected_curve_id];
    std::stringstream sstr1, sstr2;
    sstr1 << std::cout << base_curve;
    sstr2 << std::cout << exp_curve;
    if (sstr1.str() != sstr2.str()) {
      std::cout << "Obtained result and expected result does not match"
                << std::endl;
      std::cout << std::endl << "Result obtained: " << sstr1.str() << std::endl;
      std::cout << std::endl << "Expected result: " << sstr2.str() << std::endl;
      return false;
    }
#endif
  }
  else if (type == 1) {
    X_monotone_curve_2 base_curve = this->m_xcurves[id1];
    X_monotone_subcurve_2 x_segment = this->m_xsegments[segment_id];

    std::cout << "Test: push_front ( "
              << x_segment << "into"
              << base_curve << " ) ? ";

    this->m_geom_traits.push_front_2_object()(base_curve, x_segment);

    X_monotone_curve_2 exp_curve = this->m_xcurves[expected_curve_id];

    if (!this->compare_curves(exp_curve, base_curve)) return false;
  }

  else {
    std::cout << "Incorrect type of operator. "
              << "Please refer to the descriptopn in the data file."
              << std::endl;
    return false;
  }

  return true;
}

/* Compare_x_2 for polycurve
 * This functor compare_x_2 in polylines/polycurves also supports segments and
 * not just points.
 * This wrapper will only test for the x-monotone segments. For testing the
 * points, compare_x_wrapper can be used.
 */
template <typename GeomTraits>
bool Traits_test<GeomTraits>::
compare_x_wrapper(std::istringstream& str_stream) {
  unsigned int id1, id2;
  unsigned int expected_answer;
  unsigned int real_answer;

  str_stream >> id1;
  std::pair<Enum_type, unsigned int> next_input =
    this->get_next_input(str_stream);
  if (next_input.first == Base::NUMBER) {
    id2 = next_input.second;
    expected_answer = this->get_expected_enum(str_stream);
    std::cout << "Test: compare_x( "
              << this->m_points[id1] << ", "
              << this->m_points[id2] << " ) ? ";
    real_answer = this->m_geom_traits.compare_x_2_object()(this->m_points[id1],
                                                           this->m_points[id2]);
  }
  else {
    assert(next_input.first == Base::CURVE_END);
    CGAL::Arr_curve_end end1 =
      static_cast<CGAL::Arr_curve_end>(next_input.second);
    str_stream >> id2;
    next_input = this->get_next_input(str_stream);
    assert(next_input.first == Base::CURVE_END);
    CGAL::Arr_curve_end end2 =
      static_cast<CGAL::Arr_curve_end>(next_input.second);
    expected_answer = this->get_expected_enum(str_stream);
    std::cout << "Test: compare_x( "
              << this->m_xsegments[id1] << ", "
              << this->curve_end_str(end1) << ", "
              << this->m_xsegments[id2] << ", "
              << this->curve_end_str(end2) << " ) ? ";
    real_answer =
      this->m_geom_traits.compare_x_2_object()(this->m_xsegments[id1], end1,
                                               this->m_xsegments[id2], end2);
  }
  std::cout <<
    ((expected_answer == static_cast<unsigned int>(CGAL::SMALLER)) ? "SMALLER" :
     ((expected_answer == static_cast<unsigned int>(CGAL::LARGER)) ? "LARGER" :
      "EQUAL")) << " ";
  return this->compare(expected_answer, real_answer);
}

template <typename GeomTraits>
bool Traits_test<GeomTraits>::
compare_xy_wrapper(std::istringstream& str_stream) {
  unsigned int id1, id2;
  unsigned int expected_answer;
  unsigned int real_answer;

  str_stream >> id1;
  std::pair<Enum_type, unsigned int> next_input =
    this->get_next_input(str_stream);
  if (next_input.first == Base::NUMBER) {
    id2 = next_input.second;
    expected_answer = this->get_expected_enum(str_stream);
    std::cout << "Test: compare_xy( "
              << this->m_points[id1] << ", "
              << this->m_points[id2] << " ) ? ";
    real_answer =
      this->m_geom_traits.compare_xy_2_object()(this->m_points[id1],
                                                this->m_points[id2]);
  }
  else {
    assert(next_input.first == Base::CURVE_END);
    CGAL::Arr_curve_end end1 =
      static_cast<CGAL::Arr_curve_end>(next_input.second);
    str_stream >> id2;
    next_input = this->get_next_input(str_stream);
    assert(next_input.first == Base::CURVE_END);
    CGAL::Arr_curve_end end2 =
      static_cast<CGAL::Arr_curve_end>(next_input.second);
    expected_answer = this->get_expected_enum(str_stream);
    std::cout << "Test: compare_xy( "
              << this->m_xsegments[id1] << ", "
              << this->curve_end_str(end1) << ", "
              << this->m_xsegments[id2] << ", "
              << this->curve_end_str(end2) << " ) ? ";
    real_answer =
      this->m_geom_traits.compare_xy_2_object()(this->m_xsegments[id1], end1,
                                                this->m_xsegments[id2], end2);
  }

  std::cout <<
    ((expected_answer == static_cast<unsigned int>(CGAL::SMALLER)) ? "SMALLER" :
     ((expected_answer == static_cast<unsigned int>(CGAL::LARGER)) ? "LARGER" :
      "EQUAL")) << " ";
  return this->compare(expected_answer, real_answer);
}

template <typename GeomTraits>
bool Traits_test<GeomTraits>::
number_of_points_wrapper(std::istringstream& str_stream) {
  using Geom_traits = GeomTraits;
  using size_type = typename Geom_traits::size_type;

  unsigned int id;
  size_type expected_result;
  str_stream >> id >> expected_result;
  std::cout << "Test: Number_of_points( " << this->m_curves[id] << " ) ? " ;
  size_type real_answer =
    this->m_geom_traits.number_of_points_2_object()(this->m_curves[id]);
  return this->compare(expected_result, real_answer);
}

template <typename GeomTraits>
bool Traits_test<GeomTraits>::
compare_endpoints_xy_wrapper(std::istringstream& str_stream) {
  unsigned int id;
  str_stream >> id;

  unsigned int expected_answer = this->get_expected_enum(str_stream);

  std::cout << "Test: compare_endpoints_xy( " << this->m_xcurves[id]
            << " ) ? " << expected_answer << " ";
  unsigned int real_answer =
    this->m_geom_traits.compare_endpoints_xy_2_object()(this->m_xcurves[id]);

  return this->compare(expected_answer, real_answer);
}

template <typename GeomTraits>
bool Traits_test<GeomTraits>::
construct_opposite_wrapper(std::istringstream& str_stream) {
  unsigned int id1, id2;
  str_stream >> id1 >> id2;

  std::cout << "Test: construct_opposite( " << this->m_xcurves[id1] << " ) ? "
            << "expected_answer: " << this->m_xcurves[id2]<< " ";

  X_monotone_curve_2 obtained_curve =
    this->m_geom_traits.construct_opposite_2_object()(this->m_xcurves[id1]);

  return this->compare_curves(obtained_curve, this->m_xcurves[id2]);
}

#endif
//  end of POLYCURVE_CONIC_GEOM_TRAITS preprocessor if

#if TEST_GEOM_TRAITS != POLYCURVE_CONIC_GEOM_TRAITS && \
    TEST_GEOM_TRAITS != POLYCURVE_CIRCULAR_ARC_GEOM_TRAITS && \
    TEST_GEOM_TRAITS != POLYCURVE_BEZIER_GEOM_TRAITS && \
    TEST_GEOM_TRAITS != POLYLINE_GEOM_TRAITS
/*! Test Compare_x_2
 */
template <typename GeomTraits>
bool Traits_test<GeomTraits>::
compare_x_wrapper(std::istringstream& str_stream) {
  unsigned int id1, id2;
  str_stream >> id1 >> id2;
  unsigned int exp_answer = this->get_expected_enum(str_stream);
  std::cout << "Test: compare_x( " << this->m_points[id1] << ", "
            << this->m_points[id2] << " ) ? " << exp_answer << " ";

  unsigned int real_answer =
    this->m_geom_traits.compare_x_2_object()(this->m_points[id1],
                                             this->m_points[id2]);
  return this->compare(exp_answer, real_answer);
}

/*! Test Compare_xy_2
 */
template <typename GeomTraits>
bool Traits_test<GeomTraits>::
compare_xy_wrapper(std::istringstream& str_stream) {
  unsigned int id1, id2;
  str_stream >> id1 >> id2;
  unsigned int exp_answer = this->get_expected_enum(str_stream);
  std::cout << "Test: compare_xy( " << this->m_points[id1] << ", "
            << this->m_points[id2] << " ) ? " << exp_answer << " ";

  unsigned int real_answer =
    this->m_geom_traits.compare_xy_2_object()(this->m_points[id1],
                                              this->m_points[id2]);
  return this->compare(exp_answer, real_answer);
}

#endif

/*! Tests Construct_min_vertex_2.
 * Degenerate case: vertical curve.
 */
template <typename GeomTraits>
bool Traits_test<GeomTraits>::
min_vertex_wrapper(std::istringstream& str_stream) {
  unsigned int id1, id2;
  str_stream >> id1 >> id2;
  Point_2& exp_answer = this->m_points[id2];
  std::cout << "Test: min_vertex( " << CGAL::IO::oformat(this->m_xcurves[id1]) << " ) ? "
            << exp_answer << " ";

  Point_2 real_answer =
    this->m_geom_traits.construct_min_vertex_2_object()(this->m_xcurves[id1]);
  return this->compare_points(exp_answer, real_answer);
}

/*! Tests Construct_max_vertex_2.
 * Degenerate case: vertical curve.
 */
template <typename GeomTraits>
bool Traits_test<GeomTraits>::
max_vertex_wrapper(std::istringstream& str_stream) {
  unsigned int id1, id2;
  str_stream >> id1 >> id2;
  Point_2& exp_answer = this->m_points[id2];
  std::cout << "Test: max_vertex( " << CGAL::IO::oformat(this->m_xcurves[id1]) << " ) ? "
            << exp_answer << " ";

  Point_2 real_answer =
    this->m_geom_traits.construct_max_vertex_2_object()(this->m_xcurves[id1]);
  return this->compare_points(exp_answer, real_answer);
}

template <typename GeomTraits>
bool
Traits_test<GeomTraits>::
is_vertical_wrapper(std::istringstream& str_stream) {
  unsigned int id;
  str_stream >> id;
  bool exp_answer = this->get_expected_boolean(str_stream);
  std::cout << "Test: is_vertical( " << CGAL::IO::oformat(this->m_xcurves[id]) << " ) ? "
            << exp_answer << " ";

  bool real_answer =
    this->m_geom_traits.is_vertical_2_object()(this->m_xcurves[id]);
  return this->compare(exp_answer, real_answer);
}

/*! Tests Compare_y_at_x_2.
 * Return the location of the given point with respect to the input curve.
 * Degenerate cases: The point is an endpoint of the curve.
 *                   The curve is vertical.
 */
template <typename GeomTraits>
bool
Traits_test<GeomTraits>::
compare_y_at_x_wrapper(std::istringstream& str_stream) {
  unsigned int id1, id2;
  str_stream >> id1 >> id2;
  unsigned int exp_answer = this->get_expected_enum(str_stream);
  std::cout << "Test: compare_y_at_x( " << this->m_points[id1] << ","
           << CGAL::IO::oformat(this->m_xcurves[id2]) << " ) ? " << exp_answer << " ";

  unsigned int real_answer =
    this->m_geom_traits.compare_y_at_x_2_object()(this->m_points[id1],
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
template <typename GeomTraits>
bool
Traits_test<GeomTraits>::
compare_y_at_x_left_wrapper(std::istringstream& str_stream) {
  using Has_left_category = typename GeomTraits::Has_left_category;
  return compare_y_at_x_left_wrapper_imp(str_stream, Has_left_category());
}

template <typename GeomTraits>
bool
Traits_test<GeomTraits>::
compare_y_at_x_left_wrapper_imp(std::istringstream&, CGAL::Tag_false) {
  CGAL_error();
  return false;
}

template <typename GeomTraits>
bool
Traits_test<GeomTraits>::
compare_y_at_x_left_wrapper_imp(std::istringstream& str_stream, CGAL::Tag_true) {
  unsigned int id1, id2, id3;
  str_stream >> id1 >> id2 >> id3;
  unsigned int exp_answer = this->get_expected_enum(str_stream);
  std::cout << "Test: compare_y_at_x_left( " << this->m_xcurves[id1] << ","
            << this->m_xcurves[id2] << ", " << this->m_points[id3] << " ) ? "
            << exp_answer << " ";

  unsigned int real_answer = this->m_geom_traits.compare_y_at_x_left_2_object()
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
template <typename GeomTraits>
bool Traits_test<GeomTraits>::
compare_y_at_x_right_wrapper(std::istringstream& str_stream) {
  unsigned int id1, id2, id3;
  str_stream >> id1 >> id2 >> id3;
  unsigned int exp_answer = this->get_expected_enum(str_stream);
  std::cout << "Test: compare_y_at_x_right( " << CGAL::IO::oformat(this->m_xcurves[id1]) << ","
            << CGAL::IO::oformat(this->m_xcurves[id2]) << ", " << this->m_points[id3] << " ) ? "
            << exp_answer << " ";

  unsigned int real_answer =
    this->m_geom_traits.compare_y_at_x_right_2_object()
    (this->m_xcurves[id1], this->m_xcurves[id2], this->m_points[id3]);
  return this->compare(exp_answer, real_answer);
}

/*! Tests Equal_2::operator()(Point_2, Point_2).
 * Check whether two points are the same.
 */
template <typename GeomTraits>
bool Traits_test<GeomTraits>::
equal_points_wrapper(std::istringstream& str_stream) {
  unsigned int id1, id2;
  str_stream >> id1 >> id2;
  bool exp_answer = this->get_expected_boolean(str_stream);
  std::cout << "Test: equal( " << this->m_points[id1] << ", "
            << this->m_points[id2] << " ) ? " << exp_answer << " ";

  bool real_answer = this->m_geom_traits.equal_2_object()(this->m_points[id1],
                                                          this->m_points[id2]);
  return this->compare(exp_answer, real_answer);
}

/*! Tests Equal_2::operator()(X_monotone_curve_2, X_monotone_curve_2).
 * Check whether two x-monotone curves are the same.
 */
template <typename GeomTraits>
bool Traits_test<GeomTraits>::
equal_curves_wrapper(std::istringstream& str_stream) {
  unsigned int id1, id2;
  str_stream >> id1 >> id2;
  bool exp_answer = this->get_expected_boolean(str_stream);
  std::cout << "Test: equal( " << CGAL::IO::oformat(this->m_xcurves[id1]) << ", "
            << CGAL::IO::oformat(this->m_xcurves[id2]) << " ) ? " << exp_answer << " ";

  bool real_answer = this->m_geom_traits.equal_2_object()(this->m_xcurves[id1],
                                                          this->m_xcurves[id2]);
  return this->compare(exp_answer, real_answer);
}

/*! Tests Make_x_monotone_2.
 * Cut the given curve into x-monotone subcurves and insert them into the
 * given output iterator.
 * Degenerate cases for polylines: The first segment is vertical. The last
 * segment is vertical. Both first and last are vertical. An internal segment
 * is vertical.
 */
template <typename GeomTraits>
bool Traits_test<GeomTraits>::
make_x_monotone_wrapper(std::istringstream& str_stream)
{
  typedef GeomTraits                                 Traits;
  typedef typename Traits::Point_2                      Point_2;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;
  typedef std::variant<Point_2, X_monotone_curve_2>   Make_x_monotone_result;
  CGAL_USE_TYPE(typename Traits::Curve_2);

  unsigned int id;
  str_stream >> id;
  std::cout << "Test: make_x_monotone( " << CGAL::IO::oformat(this->m_curves[id]) << " ) ? ";
  std::vector<Make_x_monotone_result> objs;
  this->m_geom_traits.make_x_monotone_2_object()(this->m_curves[id],
                                                 std::back_inserter(objs));

  size_t num;
  str_stream >> num;
  if (!this->compare(num, objs.size(), "size")) return false;

  for (size_t i = 0; i < num; ++i) {
    unsigned int type;                  // 0 - point, 1 - x-monotone curve
    str_stream >> type;

    unsigned int id;                    // The id of the point or x-monotone
    str_stream >> id;                   // ... curve respectively

    const auto* xcv_ptr = std::get_if<X_monotone_curve_2>(&(objs[i]));
    if (xcv_ptr != nullptr) {
      if (!this->compare(type, 1u, "type")) return false;
      if (!this->compare_curves(this->m_xcurves[id], *xcv_ptr)) return false;
      continue;
    }

    const auto* pt_ptr = std::get_if<Point_2>(&(objs[i]));
    assert(pt_ptr != nullptr);
    if (!this->compare(type, 0u, "type")) return false;
    if (!this->compare_points(this->m_points[id], *pt_ptr)) return false;
  }
  objs.clear();
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
template <typename GeomTraits>
bool Traits_test<GeomTraits>::
intersect_wrapper(std::istringstream& str_stream)
{
  typedef GeomTraits                            Traits;
  typedef typename Traits::Point_2              Point_2;
  typedef typename Traits::X_monotone_curve_2   X_monotone_curve_2;
  typedef typename Traits::Multiplicity         Multiplicity;

  typedef std::pair<Point_2, Multiplicity>      Intersection_point;
  typedef std::variant<Intersection_point, X_monotone_curve_2>
                                                Intersection_result;
  unsigned int id1, id2;
  str_stream >> id1 >> id2;
  std::vector<Intersection_result> xections;
  this->m_geom_traits.intersect_2_object()(this->m_xcurves[id1],
                                           this->m_xcurves[id2],
                                           std::back_inserter(xections));

  std::cout << "Test: intersect( " << CGAL::IO::oformat(this->m_xcurves[id1]) << ","
            << CGAL::IO::oformat(this->m_xcurves[id2]) << " ) ? ";
  size_t num;
  str_stream >> num;
  if (! this->compare(num, xections.size(), "size")) return false;

  for (size_t i = 0; i < num; ++i) {
    unsigned int type;                  // 0 - point, 1 - x-monotone curve
    str_stream >> type;
    unsigned int id;                    // The id of the point or x-monotone
    str_stream >> id;                   // ... curve respectively
    Multiplicity multiplicity;
    if (type == 0) str_stream >> multiplicity;

    unsigned int exp_type = 1;
    const X_monotone_curve_2* cv_p =
      std::get_if<X_monotone_curve_2>(&(xections[i]));

    if (cv_p != nullptr) {
      if (! this->compare(type, exp_type, "type")) return false;
      if (! this->compare_curves(this->m_xcurves[id], *cv_p)) return false;
      continue;
    }

    exp_type = 0;
    const Intersection_point* p_p =
      std::get_if<Intersection_point>(&(xections[i]));
    assert(p_p != nullptr);
    if (! this->compare(type, exp_type, "type")) return false;
    if (! this->compare_points(this->m_points[id], p_p->first)) return false;
    if (! this->compare(multiplicity, p_p->second, "multiplicity")) return false;
  }

  xections.clear();

  return true;
}

/*! Tests Split_2.
 * Split a given x-monotone curve at a given point into two sub-curves.
 * Degenerate cases for polylines: the point and a polyline internal point
 * coincides.
 */
template <typename GeomTraits>
bool Traits_test<GeomTraits>::split_wrapper(std::istringstream& str_stream) {
  typedef GeomTraits                                 Traits;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;

  unsigned int id1, id2, id3, id4;
  str_stream >> id1 >> id2 >> id3 >> id4;
  X_monotone_curve_2 cv1, cv2;
  std::cout << "Test: split( " << CGAL::IO::oformat(this->m_xcurves[id1]) << ","
            << this->m_points[id2] << " ) ? ";

  this->m_geom_traits.split_2_object()(this->m_xcurves[id1],
                                     this->m_points[id2], cv1, cv2);
  return this->compare_curves(this->m_xcurves[id3], cv1) &&
    this->compare_curves(this->m_xcurves[id4], cv2);
}

/*! Tests Are_mergeable_2.
 * Check whether it is possible to merge two given x-monotone curves.
 */
template <typename GeomTraits>
bool Traits_test<GeomTraits>::
are_mergeable_wrapper(std::istringstream& str_stream) {
  using Has_merge_category = typename GeomTraits::Has_merge_category;
  return are_mergeable_wrapper_imp(str_stream, Has_merge_category());
}

template <typename GeomTraits>
bool
Traits_test<GeomTraits>::
are_mergeable_wrapper_imp(std::istringstream&, CGAL::Tag_false) {
  std::cout << "I am at the wrong place" << std::endl;
  CGAL_error();
  return false;
}

template <typename GeomTraits>
bool Traits_test<GeomTraits>::
are_mergeable_wrapper_imp(std::istringstream& str_stream, CGAL::Tag_true) {
  unsigned int id1, id2;
  str_stream >> id1 >> id2;
  bool exp_answer = this->get_expected_boolean(str_stream);
  std::cout << "Test: are_mergeable( " << CGAL::IO::oformat(this->m_xcurves[id1]) << ", "
            << CGAL::IO::oformat(this->m_xcurves[id2]) << " ) ? " << exp_answer << " ";

  bool real_answer =
    this->m_geom_traits.are_mergeable_2_object()(this->m_xcurves[id1],
                                                 this->m_xcurves[id2]);
  return this->compare(exp_answer, real_answer);
}

/*! Tests Merge_2.
 * Merge two given x-monotone curves into a single curve.
 */
template <typename GeomTraits>
bool Traits_test<GeomTraits>::merge_wrapper(std::istringstream& str_stream) {
  typedef typename GeomTraits::Has_merge_category         Has_merge_category;
  return merge_wrapper_imp(str_stream, Has_merge_category());
}

template <typename GeomTraits>
bool Traits_test<GeomTraits>::
merge_wrapper_imp(std::istringstream&, CGAL::Tag_false) {
  CGAL_error();
  return false;
}

template <typename GeomTraits>
bool Traits_test<GeomTraits>::
merge_wrapper_imp(std::istringstream& str_stream, CGAL::Tag_true) {
  using Traits = GeomTraits                         ;
  using X_monotone_curve_2 = typename Traits::X_monotone_curve_2;

  unsigned int id1, id2, id;
  str_stream >> id1 >> id2 >> id;
  X_monotone_curve_2 cv;
  std::cout << "Test: merge( " << CGAL::IO::oformat(this->m_xcurves[id1]) << ", "
            << CGAL::IO::oformat(this->m_xcurves[id2]) << " ) ? " << CGAL::IO::oformat(this->m_xcurves[id]) << " ";

  this->m_geom_traits.merge_2_object()(this->m_xcurves[id1],
                                       this->m_xcurves[id2], cv);
  return this->compare_curves(this->m_xcurves[id], cv);
}

//! Fall threw instance
template <typename GeomTraits>
template <typename GT>
bool Traits_test<GeomTraits>::
approximate_wrapper_impl2(std::istringstream&, long) { return false; }

//! Specialized instance
template <typename GeomTraits>
template <typename GT, typename T>
bool Traits_test<GeomTraits>::
approximate_wrapper_impl2(std::istringstream& is, int) {
  std::size_t id;
  typename GeomTraits::Approximate_number_type error;
  bool l2r;
  is >> id >> error >> l2r;
  const auto& xcv = this->m_xcurves[id];
  using Ap2 = typename GeomTraits::Approximate_point_2;
  std::vector<Ap2> apoints;

  auto approx = this->m_geom_traits.approximate_2_object();
  approx(xcv, error, std::back_inserter(apoints), l2r);

  std::size_t exp_num_apoints;
  is >> exp_num_apoints;
  if (apoints.size() != exp_num_apoints) {
    std::cout << "xcv: " << xcv << std::endl;
    std::cout << "error: " << error << std::endl;
    for (const auto& pt : apoints) {
      std::cout << pt << std::endl;
    }
    std::cerr << "Error: no. of inexact points does not match ("
              << apoints.size() << ", " << exp_num_apoints
              << ")!\n";
    return false;
  }
  return true;
}

//! Fall threw instance
template <typename GeomTraits>
template <typename GT>
bool Traits_test<GeomTraits>::
approximate_wrapper_impl1(std::istringstream&, long) { return false; }

//! Specialized instance
template <typename GeomTraits>
template <typename GT, typename T>
bool Traits_test<GeomTraits>::
approximate_wrapper_impl1(std::istringstream& is, int)
{ return approximate_wrapper_impl2<GeomTraits>(is, 0); }

/*! Tests Approximate_2.
 * Return an approximation of a point coordinate.
 */
template <typename GeomTraits>
bool Traits_test<GeomTraits>::approximate_wrapper(std::istringstream& is)
{ return approximate_wrapper_impl1<GeomTraits>(is, 0); }

/*! tests Construct_x_monotone_curve_2.
 * Return an x-monotone curve connecting the two given endpoints.
 */
template <typename GeomTraits>
bool Traits_test<GeomTraits>::
construct_x_monotone_curve_wrapper(std::istringstream&) { return false; }

// ///////////////////////////////////////////////////////////////////////////
// boundary-specific functors

// ---------------------------------------------------------------------------
// left-right

/*! Test Parameter_space_in_x_2
*/
template <typename GeomTraits>
bool Traits_test<GeomTraits>::
parameter_space_in_x_wrapper(std::istringstream& str_stream) {
  using Left_side_category = typename
    CGAL::internal::Arr_complete_left_side_category<GeomTraits>::Category;
  using Right_side_category = typename
    CGAL::internal::Arr_complete_right_side_category<GeomTraits>::Category;
  using LR = CGAL::internal::Arr_left_right_implementation_dispatch
    <Left_side_category, Right_side_category>;
  using Psx_tag = typename LR::Parameter_space_in_x_2_curve_end_tag;
  return parameter_space_in_x_wrapper_imp(str_stream, Psx_tag());
}

template <typename GeomTraits>
bool Traits_test<GeomTraits>::
parameter_space_in_x_wrapper_imp(std::istringstream&, CGAL::Arr_use_dummy_tag) {
  CGAL_error();
  return false;
}

template <typename GeomTraits>
bool Traits_test<GeomTraits>::
parameter_space_in_x_wrapper_imp(std::istringstream& str_stream,
                                 CGAL::Arr_use_traits_tag) {
  CGAL::Arr_parameter_space exp_answer, real_answer;
  unsigned int id;
  str_stream >> id;
  std::pair<Enum_type, unsigned int> next_input =
    this->get_next_input(str_stream);
  bool curves_op = (next_input.first == Base::CURVE_END);
  CGAL::Arr_curve_end cv_end = CGAL::ARR_MIN_END;
  if (curves_op) {
    cv_end = static_cast<CGAL::Arr_curve_end>(next_input.second);
    std::cout << "Test: parameter_space_x( " << CGAL::IO::oformat(this->m_xcurves[id]) << " , "
              << (cv_end == CGAL::ARR_MIN_END ? "MIN_END" : "MAX_END")
              << " ) ? ";
    next_input = this->get_next_input(str_stream);
    assert(next_input.first == Base::PARAMETER_SPACE);
  }
  else if (next_input.first == Base::PARAMETER_SPACE) {
    std::cout << "Test: parameter_space_in_x( " << this->m_points[id]
              << " ) ? ";
  }
  else {
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
    this->m_geom_traits.parameter_space_in_x_2_object()(this->m_xcurves[id],
                                                        cv_end) :
    this->m_geom_traits.parameter_space_in_x_2_object()(this->m_points[id]);
  return this->compare(exp_answer, real_answer);
}

/*! Test Compare_y_near_boundary_2
 */
template <typename GeomTraits>
bool Traits_test<GeomTraits>::
compare_y_near_boundary_wrapper(std::istringstream& str_stream) {
  using Left_side_category = typename
    CGAL::internal::Arr_complete_left_side_category<GeomTraits >::Category;
  using Right_side_category = typename
    CGAL::internal::Arr_complete_right_side_category<GeomTraits >::Category;
  using LR = CGAL::internal::Arr_left_right_implementation_dispatch
    <Left_side_category, Right_side_category>;
  using Cmp_tag = typename LR::Compare_y_near_boundary_2_curve_ends_tag;
  return compare_y_near_boundary_wrapper_imp(str_stream, Cmp_tag());
}

template <typename GeomTraits>
bool Traits_test<GeomTraits>::
compare_y_near_boundary_wrapper_imp(std::istringstream&,
                                    CGAL::Arr_use_dummy_tag) {
  CGAL_error();
  return false;
}

template <typename GeomTraits>
bool Traits_test<GeomTraits>::
compare_y_near_boundary_wrapper_imp(std::istringstream& str_stream,
                                    CGAL::Arr_use_traits_tag) {
  unsigned int id1, id2;
  str_stream >> id1 >> id2;
  std::pair<Enum_type, unsigned int> next_input =
    this->get_next_input(str_stream);
  assert(next_input.first == Base::CURVE_END);
  CGAL::Arr_curve_end cv_end =
    static_cast<CGAL::Arr_curve_end>(next_input.second);
  std::cout << "Test: compare_y_near_boundary( " << CGAL::IO::oformat(this->m_xcurves[id1])
            << " , " << CGAL::IO::oformat(this->m_xcurves[id2]) << " , "
            << (cv_end == CGAL::ARR_MIN_END ? "MIN_END" : "MAX_END")
            << " ) ? ";
  next_input = this->get_next_input(str_stream);
  assert(next_input.first == Base::SIGN);
  CGAL::Comparison_result exp_answer =
    static_cast<CGAL::Comparison_result>(next_input.second);
  std::cout << (exp_answer == CGAL::SMALLER ? "SMALLER":
               (exp_answer == CGAL::LARGER ? "LARGER" : "EQUAL")) << " ";

  CGAL::Comparison_result real_answer =
    this->m_geom_traits.compare_y_near_boundary_2_object()
    (this->m_xcurves[id1], this->m_xcurves[id2], cv_end);
  return this->compare(exp_answer, real_answer);
}

// TODO Is_on_y_identification_2

// ---------------------------------------------------------------------------
// bottom-top

/*! Test Parameter_space_in_y_2
 */
template <typename GeomTraits>
bool Traits_test<GeomTraits>::
parameter_space_in_y_wrapper(std::istringstream& str_stream) {
  using Bottom_side_category = typename
    CGAL::internal::Arr_complete_bottom_side_category<GeomTraits>::Category;
  using Top_side_category = typename
    CGAL::internal::Arr_complete_top_side_category<GeomTraits>::Category;
  using BT = CGAL::internal::Arr_bottom_top_implementation_dispatch
    <Bottom_side_category, Top_side_category>;
  using Psy_tag = typename BT::Parameter_space_in_y_2_curve_end_tag;
  return parameter_space_in_y_wrapper_imp(str_stream, Psy_tag());
}

template <typename GeomTraits>
bool Traits_test<GeomTraits>::
parameter_space_in_y_wrapper_imp(std::istringstream&, CGAL::Arr_use_dummy_tag) {
  CGAL_error();
  return false;
}

template <typename GeomTraits>
bool Traits_test<GeomTraits>::
parameter_space_in_y_wrapper_imp(std::istringstream& str_stream,
                                 CGAL::Arr_use_traits_tag) {
  unsigned int id;
  str_stream >> id;
  std::pair<Enum_type, unsigned int> next_input =
    this->get_next_input(str_stream);
  bool curves_op = (next_input.first == Base::CURVE_END);
  CGAL::Arr_curve_end cv_end = CGAL::ARR_MIN_END;
  if (curves_op) {
    cv_end = static_cast<CGAL::Arr_curve_end>(next_input.second);
    std::cout << "Test: parameter_space_y( " << CGAL::IO::oformat(this->m_xcurves[id]) << " , "
              << (cv_end == CGAL::ARR_MIN_END ? "MIN_END" : "MAX_END")
              << " ) ? ";
    next_input = this->get_next_input(str_stream);
    assert(next_input.first == Base::PARAMETER_SPACE);
  }
  else if (next_input.first == Base::PARAMETER_SPACE) {
    std::cout << "Test: parameter_space_in_y( " << this->m_points[id]
              << " ) ? ";
  }
  else {
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
    this->m_geom_traits.parameter_space_in_y_2_object()(this->m_xcurves[id],
                                                        cv_end) :
    this->m_geom_traits.parameter_space_in_y_2_object()(this->m_points[id]);
  return this->compare(exp_answer, real_answer);
}

/* Compare_x_near_boundary_2
 */
template <typename GeomTraits>
bool Traits_test<GeomTraits>::
compare_x_near_boundary_wrapper(std::istringstream& str_stream) {
  using Bottom_side_category = typename
    CGAL::internal::Arr_complete_bottom_side_category<GeomTraits>::Category;
  using Top_side_category = typename
    CGAL::internal::Arr_complete_top_side_category<GeomTraits>::Category;
  using BT = CGAL::internal::Arr_bottom_top_implementation_dispatch
    <Bottom_side_category, Top_side_category>;
  using Cmp_tag = typename BT::Compare_x_near_boundary_2_curve_ends_tag;
  return compare_x_near_boundary_wrapper_imp(str_stream, Cmp_tag());
}

template <typename GeomTraits>
bool Traits_test<GeomTraits>::
compare_x_near_boundary_wrapper_imp(std::istringstream&,
                                    CGAL::Arr_use_dummy_tag) {
  CGAL_error();
  return false;
}

template <typename GeomTraits>
bool Traits_test<GeomTraits>::
compare_x_near_boundary_wrapper_imp(std::istringstream& str_stream,
                                    CGAL::Arr_use_traits_tag) {
  std::cout << "Test: compare_x_near_boundary( ";

  unsigned int id1, id2;
  // xcurve xcurve curve-end
  str_stream >> id1 >> id2;
  std::pair<Enum_type, unsigned int> next_input =
    this->get_next_input(str_stream);
  assert(next_input.first == Base::CURVE_END);
  CGAL::Arr_curve_end cv_end =
    static_cast<CGAL::Arr_curve_end>(next_input.second);

  std::cout << CGAL::IO::oformat(this->m_xcurves[id1]) << ", "
            << CGAL::IO::oformat(this->m_xcurves[id2]) << " , "
            << this->curve_end_str(cv_end) << " ) ? ";

  next_input = this->get_next_input(str_stream);
  assert(next_input.first == Base::SIGN);
  CGAL::Comparison_result exp_answer =
    static_cast<CGAL::Comparison_result>(next_input.second);
  std::cout << (exp_answer == CGAL::SMALLER ? "SMALLER":
               (exp_answer == CGAL::LARGER ? "LARGER":"EQUAL")) << " ";

  CGAL::Comparison_result real_answer =
    this->m_geom_traits.compare_x_near_boundary_2_object()
    (this->m_xcurves[id1], this->m_xcurves[id2], cv_end);
  return this->compare(exp_answer, real_answer);
}

/* Compare_x_on_boundary
 */
template <typename GeomTraits>
bool Traits_test<GeomTraits>::
compare_x_on_boundary_wrapper(std::istringstream& str_stream) {
  using Bottom_side_category = typename
    CGAL::internal::Arr_complete_bottom_side_category<GeomTraits>::Category;
  using Top_side_category = typename
    CGAL::internal::Arr_complete_top_side_category<GeomTraits>::Category;
  using BT = CGAL::internal::Arr_bottom_top_implementation_dispatch
    <Bottom_side_category, Top_side_category>;
  using Cmp_tag1 = typename BT::Compare_x_on_boundary_2_points_tag;
  using Cmp_tag2 = typename BT::Compare_x_on_boundary_2_point_curve_end_tag;
  using Cmp_tag3 = typename BT::Compare_x_on_boundary_2_curve_ends_tag;
  using Cmp_tag12 = typename CGAL::internal::Or_traits<Cmp_tag1, Cmp_tag2>::type;
  using Cmp_tag = typename CGAL::internal::Or_traits<Cmp_tag12, Cmp_tag3>::type;
  return compare_x_on_boundary_wrapper_imp(str_stream, Cmp_tag());
}

template <typename GeomTraits>
bool Traits_test<GeomTraits>::
compare_x_on_boundary_wrapper_imp(std::istringstream&, CGAL::Arr_use_dummy_tag) {
  CGAL_error();
  return false;
}

template <typename GeomTraits>
bool Traits_test<GeomTraits>::
compare_x_on_boundary_wrapper_imp(std::istringstream& str_stream,
                                  CGAL::Arr_use_traits_tag) {
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
    assert(next_input.first == Base::CURVE_END);
    cv_end1 = static_cast<CGAL::Arr_curve_end>(next_input.second);

    std::cout << this->m_points[id1] << " , "
              << this->m_xcurves[id2]<< " , "
              << this->curve_end_str(cv_end1) << " ) ? ";
  }
  else if (next_input.first == Base::CURVE_END) {
    cv_end1 = static_cast<CGAL::Arr_curve_end>(next_input.second);
    next_input = this->get_next_input(str_stream);
    assert(next_input.first == Base::NUMBER);
    id2 = static_cast<unsigned int>(next_input.second);
    next_input = this->get_next_input(str_stream);
    assert(next_input.first == Base::CURVE_END);
    cv_end2 = static_cast<CGAL::Arr_curve_end>(next_input.second);

    std::cout << this->m_xcurves[id1] << ", "
              << this->curve_end_str(cv_end1) << ", "
              << this->m_xcurves[id2] << " , "
              << this->curve_end_str(cv_end2) << " ) ? ";
  }
  else CGAL_error();

  next_input = this->get_next_input(str_stream);
  assert(next_input.first == Base::SIGN);
  CGAL::Comparison_result exp_answer =
    static_cast<CGAL::Comparison_result>(next_input.second);
  std::cout << (exp_answer == CGAL::SMALLER ? "SMALLER":
               (exp_answer == CGAL::LARGER ? "LARGER":"EQUAL")) << " ";

  real_answer = (curves_op) ?
    this->m_geom_traits.compare_x_on_boundary_2_object()
    (this->m_points[id1], this->m_xcurves[id2], cv_end1) :
    this->m_geom_traits.compare_x_on_boundary_2_object()
    (this->m_xcurves[id1], cv_end1,this->m_xcurves[id2], cv_end2);
  return this->compare(exp_answer, real_answer);
}

// TODO Is_on_x_identification_2

#endif
