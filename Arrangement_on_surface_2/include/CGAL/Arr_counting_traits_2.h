// Copyright (c) 2005,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Efi Fogel      <efif@post.tau.ac.il>
//            Eric Berberich <ericb@post.tau.ac.il>

#ifndef CGAL_ARR_COUNTING_TRAITS_H
#define CGAL_ARR_COUNTING_TRAITS_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>

/*! \file
 * A counting traits-class for the arrangement package.
 * This is a meta-traits class. It is parameterized with another traits class
 * and inherits from it. For each traits method it maintains a counter that
 * counts the number of invokations into the method.
 */

#include <iostream>
#include <string.h>
#include <atomic>
#include <array>

#include <CGAL/basic.h>
#include <CGAL/Arr_enums.h>
#include <CGAL/Arr_tags.h>

namespace CGAL {

/*! \class
 * A model of the ArrangementTraits_2 concept that counts the methods invoked.
 */
template <typename Base_traits>
class Arr_counting_traits_2 : public Base_traits {
public:
  enum Operation_id {
    COMPARE_X_OP = 0,
    COMPARE_XY_OP,
    CONSTRUCT_MIN_VERTEX_OP,
    CONSTRUCT_MAX_VERTEX_OP,
    IS_VERTICAL_OP,
    COMPARE_Y_AT_X_OP,
    EQUAL_POINTS_OP,
    EQUAL_CURVES_OP,
    COMPARE_Y_AT_X_LEFT_OP,
    COMPARE_Y_AT_X_RIGHT_OP,
    MAKE_X_MONOTONE_OP,
    SPLIT_OP,
    INTERSECT_OP,
    ARE_MERGEABLE_OP,
    MERGE_OP,
    CONSTRUCT_OPPOSITE_OP,
    COMPARE_ENDPOINTS_XY_OP,

    PARAMETER_SPACE_IN_X_CURVE_END_OP,
    PARAMETER_SPACE_IN_X_POINT_OP,
    PARAMETER_SPACE_IN_X_CURVE_OP,
    IS_ON_X_IDENTIFICATION_POINT_OP,
    IS_ON_X_IDENTIFICATION_CURVE_OP,
    COMPARE_Y_ON_BOUNDARY_OP,
    COMPARE_Y_NEAR_BOUNDARY_OP,

    PARAMETER_SPACE_IN_Y_CURVE_END_OP,
    PARAMETER_SPACE_IN_Y_POINT_OP,
    PARAMETER_SPACE_IN_Y_CURVE_OP,
    IS_ON_Y_IDENTIFICATION_POINT_OP,
    IS_ON_Y_IDENTIFICATION_CURVE_OP,
    COMPARE_X_ON_BOUNDARY_POINTS_OP,
    COMPARE_X_ON_BOUNDARY_POINT_CURVE_END_OP,
    COMPARE_X_ON_BOUNDARY_CURVE_ENDS_OP,
    COMPARE_X_NEAR_BOUNDARY_OP,

    NUMBER_OF_OPERATIONS
  };

  typedef Base_traits                           Base;
  typedef Arr_counting_traits_2<Base>           Self;

  /*! Construct default */
  Arr_counting_traits_2() : Base()
  {
    clear_counters();
    increment();
  }

  /*! Construct copy */
  Arr_counting_traits_2(const Arr_counting_traits_2& other) : Base(other)
  {
    clear_counters();
    increment();
  }

  /*! Obtain the counter of the given operation */
  size_t count(Operation_id id) const
  { return m_counters[id]; }

  size_t count_compare_x() const
  { return m_counters[COMPARE_X_OP]; }

  size_t count_compare_xy() const
  { return m_counters[COMPARE_XY_OP]; }

  size_t count_construct_min_vertex() const
  { return m_counters[CONSTRUCT_MIN_VERTEX_OP]; }

  size_t count_construct_max_vertex() const
  { return m_counters[CONSTRUCT_MAX_VERTEX_OP]; }

  size_t count_is_vertical() const
  { return m_counters[IS_VERTICAL_OP]; }

  size_t count_compare_y_at_x() const
  { return m_counters[COMPARE_Y_AT_X_OP]; }

  size_t count_equal_points() const
  { return m_counters[EQUAL_POINTS_OP]; }

  size_t count_equal_curves() const
  { return m_counters[EQUAL_CURVES_OP]; }

  size_t count_compare_y_at_x_left() const
  { return m_counters[COMPARE_Y_AT_X_LEFT_OP]; }

  size_t count_compare_y_at_x_right() const
  { return m_counters[COMPARE_Y_AT_X_RIGHT_OP]; }

  size_t count_make_x_monotone() const
  { return m_counters[MAKE_X_MONOTONE_OP]; }

  size_t count_split() const
  { return m_counters[SPLIT_OP]; }

  size_t count_intersect() const
  { return m_counters[INTERSECT_OP]; }

  size_t count_are_mergeable() const
  { return m_counters[ARE_MERGEABLE_OP]; }

  size_t count_merge() const
  { return m_counters[MERGE_OP]; }

  size_t count_construct_opposite() const
  { return m_counters[CONSTRUCT_OPPOSITE_OP]; }

  size_t count_compare_endpoints_xy() const
  { return m_counters[COMPARE_ENDPOINTS_XY_OP]; }

  // left-right

  size_t count_parameter_space_in_x_curve_end() const
  { return m_counters[PARAMETER_SPACE_IN_X_CURVE_END_OP]; }

  size_t count_parameter_space_in_x_curve() const
  { return m_counters[PARAMETER_SPACE_IN_X_CURVE_OP]; }

  size_t count_parameter_space_in_x_point() const
  { return m_counters[PARAMETER_SPACE_IN_X_POINT_OP]; }

  size_t count_is_on_x_identification_point() const
  { return m_counters[IS_ON_X_IDENTIFICATION_POINT_OP]; }

  size_t count_is_on_x_identification_curve() const
  { return m_counters[IS_ON_X_IDENTIFICATION_CURVE_OP]; }

  size_t count_compare_y_on_boundary() const
  { return m_counters[COMPARE_Y_ON_BOUNDARY_OP]; }

  size_t count_compare_y_near_boundary() const
  { return m_counters[COMPARE_Y_NEAR_BOUNDARY_OP]; }


  // bottom-top

  size_t count_parameter_space_in_y_curve_end() const
  { return m_counters[PARAMETER_SPACE_IN_Y_CURVE_END_OP]; }

  size_t count_parameter_space_in_y_curve() const
  { return m_counters[PARAMETER_SPACE_IN_Y_CURVE_OP]; }

  size_t count_parameter_space_in_y_point() const
  { return m_counters[PARAMETER_SPACE_IN_Y_POINT_OP]; }

  size_t count_is_on_y_identification_point() const
  { return m_counters[IS_ON_Y_IDENTIFICATION_POINT_OP]; }

  size_t count_is_on_y_identification_curve() const
  { return m_counters[IS_ON_Y_IDENTIFICATION_CURVE_OP]; }

  size_t count_compare_x_on_boundary_points() const
  { return m_counters[COMPARE_X_ON_BOUNDARY_POINTS_OP]; }

  size_t count_compare_x_on_boundary_point_curve_end() const
  { return m_counters[COMPARE_X_ON_BOUNDARY_POINT_CURVE_END_OP]; }

  size_t count_compare_x_on_boundary_curve_ends() const
  { return m_counters[COMPARE_X_ON_BOUNDARY_CURVE_ENDS_OP]; }

  size_t count_compare_x_near_boundary() const
  { return m_counters[COMPARE_X_NEAR_BOUNDARY_OP]; }

  /// \name Types and functors inherited from the base
  //@{

  // Traits types:
  typedef typename Base::Has_left_category          Has_left_category;
  typedef typename Base::Has_merge_category         Has_merge_category;
  typedef typename Base::Has_do_intersect_category  Has_do_intersect_category;

  typedef typename internal::Arr_complete_left_side_category< Base >::Category
                                                    Left_side_category;
  typedef typename internal::Arr_complete_bottom_side_category< Base >::Category
                                                    Bottom_side_category;
  typedef typename internal::Arr_complete_top_side_category< Base >::Category
                                                    Top_side_category;
  typedef typename internal::Arr_complete_right_side_category< Base >::Category
                                                    Right_side_category;

  typedef typename Base::Point_2                    Point_2;
  typedef typename Base::X_monotone_curve_2         X_monotone_curve_2;
  typedef typename Base::Curve_2                    Curve_2;

  /*! A functor that compares the x-coordinates of two points */
  class Compare_x_2 {
  private:
    typename Base::Compare_x_2 m_object;
    size_t& m_counter;

  public:
    /*! Construct */
    Compare_x_2(const Base* base, size_t& counter) :
      m_object(base->compare_x_2_object()), m_counter(counter) {}

    /*! Operate */
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const
    { ++m_counter; return m_object(p1, p2); }
  };

  /*! A functor that compares two points lexigoraphically: by x, then by y. */
  class Compare_xy_2 {
  private:
    typename Base::Compare_xy_2 m_object;
    size_t& m_counter;

  public:
    /*! Construct */
    Compare_xy_2(const Base* base, size_t& counter) :
      m_object(base->compare_xy_2_object()), m_counter(counter) {}

    /*! Operate */
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const
    { ++m_counter; return m_object(p1, p2); }
  };

  /*! A functor that obtains the left endpoint of an x-monotone curve. */
  class Construct_min_vertex_2 {
  private:
    typename Base::Construct_min_vertex_2 m_object;
    size_t& m_counter;

  public:
    /*! Construct */
    Construct_min_vertex_2(const Base* base, size_t& counter) :
      m_object(base->construct_min_vertex_2_object()), m_counter(counter) {}

    /*! Operate */
    const Point_2 operator()(const X_monotone_curve_2& xc) const
    { ++m_counter; return m_object(xc); }
  };

  /*! A functor that obtains the right endpoint of an x-monotone curve. */
  class Construct_max_vertex_2 {
  private:
    typename Base::Construct_max_vertex_2 m_object;
    size_t& m_counter;

  public:
    /*! Construct */
    Construct_max_vertex_2(const Base* base, size_t& counter) :
      m_object(base->construct_max_vertex_2_object()), m_counter(counter) {}

    /*! Operate */
    const Point_2 operator()(const X_monotone_curve_2& xc) const
    { ++m_counter; return m_object(xc); }
  };

  /*! A functor that checks whether a given x-monotone curve is vertical. */
  class Is_vertical_2 {
  private:
    typename Base::Is_vertical_2 m_object;
    size_t& m_counter;

  public:
    /*! Construct */
    Is_vertical_2(const Base* base, size_t& counter) :
      m_object(base->is_vertical_2_object()), m_counter(counter) {}

    /*! Operate */
    bool operator()(const X_monotone_curve_2& xc) const
    { ++m_counter; return m_object(xc); }
  };

  /*! A functor that compares the y-coordinates of a point and an
   * x-monotone curve at the point x-coordinate.
   */
  class Compare_y_at_x_2 {
  private:
    typename Base::Compare_y_at_x_2 m_object;
    size_t& m_counter;

  public:
    /*! Construct */
    Compare_y_at_x_2(const Base* base, size_t& counter) :
      m_object(base->compare_y_at_x_2_object()), m_counter(counter) {}

    /*! Operate */
    Comparison_result operator()(const Point_2& p,
                                 const X_monotone_curve_2& xc) const
    { ++m_counter; return m_object(p, xc); }
  };

  /*! A functor that checks whether two points and two x-monotone curves are
   * identical.
   */
  class Equal_2 {
  private:
    typename Base::Equal_2 m_object;
    size_t& m_counter1;
    size_t& m_counter2;

  public:
    /*! Construct */
    Equal_2(const Base* base, size_t& counter1, size_t& counter2) :
      m_object(base->equal_2_object()),
      m_counter1(counter1), m_counter2(counter2)
    {}

    /*! Operate */
    bool operator()(const X_monotone_curve_2& xc1,
                    const X_monotone_curve_2& xc2) const
    { ++m_counter1; return m_object(xc1, xc2); }

    /*! Operate */
    bool operator()(const Point_2& p1, const Point_2& p2) const
    { ++m_counter2; return m_object(p1, p2); }
  };

  /*! A functor that compares compares the y-coordinates of two x-monotone
   * curves immediately to the left of their intersection point.
   */
  class Compare_y_at_x_left_2 {
  private:
    typename Base::Compare_y_at_x_left_2 m_object;
    size_t& m_counter;

  public:
    /*! Construct */
    Compare_y_at_x_left_2(const Base* base, size_t& counter) :
      m_object(base->compare_y_at_x_left_2_object()), m_counter(counter) {}

    /*! Operate */
    Comparison_result operator()(const X_monotone_curve_2& xc1,
                                 const X_monotone_curve_2& xc2,
                                 const Point_2& p) const
    { ++m_counter; return m_object(xc1, xc2, p); }
  };

  /*! A functor that compares compares the y-coordinates of two x-monotone
   * curves immediately to the right of their intersection point.
   */
  class Compare_y_at_x_right_2 {
  private:
    typename Base::Compare_y_at_x_right_2 m_object;
    size_t& m_counter;

  public:
    /*! Construct */
    Compare_y_at_x_right_2(const Base* base, size_t& counter) :
      m_object(base->compare_y_at_x_right_2_object()), m_counter(counter) {}

    /*! Operate */
    Comparison_result operator()(const X_monotone_curve_2& xc1,
                                 const X_monotone_curve_2& xc2,
                                 const Point_2& p) const
    { ++m_counter; return m_object(xc1, xc2, p); }
  };

  /*! \class Make_x_monotone_2
   * A functor that subdivides a curve into x-monotone curves.
   */
  class Make_x_monotone_2 {
  private:
    typename Base::Make_x_monotone_2 m_object;
    size_t& m_counter;

  public:
    /*! Construct */
    Make_x_monotone_2(const Base* base, size_t& counter) :
      m_object(base->make_x_monotone_2_object()), m_counter(counter) {}

    /*! Subdivide a given curve into x-monotone subcurves and insert them into
     * a given output iterator.
     * \param cv the curve.
     * \param oi the output iterator for the result. Its value type is a variant
     *           that wraps Point_2 or an X_monotone_curve_2 objects.
     * \return The past-the-end iterator.
     */
    template <typename OutputIterator>
    OutputIterator operator()(const Curve_2& cv, OutputIterator oi) const
    { ++m_counter; return m_object(cv, oi); }
  };

  /*! A functor that splits an arc at a point. */
  class Split_2 {
  private:
    typename Base::Split_2 m_object;
    size_t& m_counter;

  public:
    /*! Construct */
    Split_2(const Base* base, size_t& counter) :
      m_object(base->split_2_object()), m_counter(counter) {}

    /*! Operate */
    void operator()(const X_monotone_curve_2& xc, const Point_2& p,
                    X_monotone_curve_2& xc1, X_monotone_curve_2& xc2) const
    { ++m_counter; m_object(xc, p, xc1, xc2); }
  };

  /*! A functor that computes intersections between x-monotone curves. */
  class Intersect_2 {
  private:
    typename Base::Intersect_2 m_object;
    size_t& m_counter;

  public:
    /*! Construct */
    Intersect_2(const Base* base, size_t& counter) :
      m_object(base->intersect_2_object()), m_counter(counter) {}

    /*! Operate */
    template <typename OutputIterator>
    OutputIterator operator()(const X_monotone_curve_2& xc1,
                              const X_monotone_curve_2& xc2,
                              OutputIterator oi) const
    { ++m_counter; return m_object(xc1, xc2, oi); }
  };

  /*! A functor that tests whether two x-monotone curves can be merged. */
  class Are_mergeable_2 {
  private:
    typename Base::Are_mergeable_2 m_object;
    size_t& m_counter;

  public:
    /*! Construct */
    Are_mergeable_2(const Base* base, size_t& counter) :
      m_object(base->are_mergeable_2_object()), m_counter(counter) {}

    /*! Operate */
    bool operator()(const X_monotone_curve_2& xc1,
                    const X_monotone_curve_2& xc2) const
    { ++m_counter; return m_object(xc1, xc2); }
  };

  /*! A functor that merges two x-monotone curves into one. */
  class Merge_2 {
  private:
    typename Base::Merge_2 m_object;
    size_t& m_counter;

  public:
    /*! Construct */
    Merge_2(const Base* base, size_t& counter) :
      m_object(base->merge_2_object()), m_counter(counter) {}

    /*! Operate */
    void operator()(const X_monotone_curve_2& xc1,
                    const X_monotone_curve_2& xc2,
                    X_monotone_curve_2& xc) const
    { ++m_counter; m_object(xc1, xc2, xc); }
  };

  /*! A fnuctor that constructs an opposite x-monotone curve. */
  class Construct_opposite_2 {
  private:
    typename Base::Construct_opposite_2 m_object;
    size_t& m_counter;

  public:
    /*! Construct */
    Construct_opposite_2(const Base* base, size_t& counter) :
      m_object(base->construct_opposite_2_object()), m_counter(counter) {}

    /*! Operate */
    X_monotone_curve_2 operator()(const X_monotone_curve_2& xc)
    { ++m_counter; return m_object(xc); }
  };

  /*! A functor that compares the two endpoints of an x-monotone curve
   * lexigoraphically.
   */
  class Compare_endpoints_xy_2 {
  private:
    typename Base::Compare_endpoints_xy_2 m_object;
    size_t& m_counter;

  public:
    /*! Construct */
    Compare_endpoints_xy_2(const Base* base, size_t& counter) :
      m_object(base->compare_endpoints_xy_2_object()), m_counter(counter) {}

    /*! Operate */
    Comparison_result operator()(const X_monotone_curve_2& xc)
    { ++m_counter; return m_object(xc); }
  };

  // left-right

  /*! A functor that determines whether an endpoint of an x-monotone curve lies
   * on a boundary of the parameter space along the x axis.
   */
  class Parameter_space_in_x_2 {
  private:
    typename Base::Parameter_space_in_x_2 m_object;
    size_t& m_counter1;
    size_t& m_counter2;
    size_t& m_counter3;

  public:
    /*! Construct */
    Parameter_space_in_x_2(const Base* base, size_t& counter1,
                           size_t& counter2, size_t& counter3) :
      m_object(base->parameter_space_in_x_2_object()),
      m_counter1(counter1),
      m_counter2(counter2),
      m_counter3(counter3)
    {}

    /*! Operate */
    Arr_parameter_space operator()(const X_monotone_curve_2& xc,
                             Arr_curve_end ce) const
    { ++m_counter1; return m_object(xc, ce); }


    /*! Operate */
    Arr_parameter_space operator()(const Point_2& p) const
    { ++m_counter2; return m_object(p); }

    /*! Operate */
    Arr_parameter_space operator()(const X_monotone_curve_2& xc) const
    { ++m_counter3; return m_object(xc); }
  };

  /*! A functor that determines whether a point or a curve lies on an
   * identification in x.
   */
  class Is_on_x_identification_2 {
  private:
    typename Base::Is_on_x_identificiation_2 m_object;
    size_t& m_counter1;
    size_t& m_counter2;

  public:
    /*! Construct */
    Is_on_x_identification_2(const Base* base,
                             size_t& counter1, size_t& counter2) :
      m_object(base->is_on_x_identificiation_2_object()),
      m_counter1(counter1),
      m_counter2(counter2)
    {}

    /*! Operate */
    Arr_parameter_space operator()(const Point_2& p) const
    { ++m_counter1; return m_object(p); }

    /*! Operate */
    Arr_parameter_space operator()(const X_monotone_curve_2& xc) const
    { ++m_counter2; return m_object(xc); }
  };

  /*! A functor that compares the y-coordinate of two given points
   * that lie on vertical boundaries.
   */
  class Compare_y_on_boundary_2 {
  private:
    typename Base::Compare_y_on_boundary_2 m_object;
    size_t& m_counter;

  public:
    /*! Construct */
    Compare_y_on_boundary_2(const Base* base, size_t& counter) :
      m_object(base->compare_y_on_boundary_2_object()),
      m_counter(counter)
    {}

    /*! Operate */
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const
    { ++m_counter; return m_object(p1, p2); }
  };

  /*! A functor that compares the y-coordinates of curve ends near the
   * boundary of the parameter space.
   */
  class Compare_y_near_boundary_2 {
  private:
    typename Base::Compare_y_near_boundary_2 m_object;
    size_t& m_counter;

  public:
    /*! Construct */
    Compare_y_near_boundary_2(const Base* base, size_t& counter) :
      m_object(base->compare_y_near_boundary_2_object()), m_counter(counter) {}

    /*! Operate */
    Comparison_result operator()(const X_monotone_curve_2& xc1,
                                 const X_monotone_curve_2& xc2,
                                 Arr_curve_end ce) const
    { ++m_counter; return m_object(xc1, xc2, ce); }
  };

  // bottom-top

  /*! A functor that determines whether an endpoint of an x-monotone arc lies
   * on a boundary of the parameter space along the y axis.
   */
  class Parameter_space_in_y_2 {
  private:
    typename Base::Parameter_space_in_y_2 m_object;
    size_t& m_counter1;
    size_t& m_counter2;
    size_t& m_counter3;

  public:
    /*! Construct */
    Parameter_space_in_y_2(const Base* base, size_t& counter1,
                           size_t& counter2, size_t& counter3) :
      m_object(base->parameter_space_in_y_2_object()),
      m_counter1(counter1),
      m_counter2(counter2),
      m_counter3(counter3)
    {}

    /*! Operate */
    Arr_parameter_space operator()(const X_monotone_curve_2& xc,
                                   Arr_curve_end ce) const
    { ++m_counter1; return m_object(xc, ce); }

    /*! Operate */
    Arr_parameter_space operator()(const Point_2& p) const
    { ++m_counter2; return m_object(p); }

    /*! Operate */
    Arr_parameter_space operator()(const X_monotone_curve_2& xc) const
    { ++m_counter3; return m_object(xc); }
  };

  /*! A functor that determines whether a point or a curve lies on an
   * identification in x.
   */
  class Is_on_y_identification_2 {
  private:
    typename Base::Is_on_y_identificiation_2 m_object;
    size_t& m_counter1;
    size_t& m_counter2;

  public:
    /*! Construct */
    Is_on_y_identification_2(const Base* base,
                             size_t& counter1, size_t& counter2) :
      m_object(base->is_on_y_identificiation_2_object()),
      m_counter1(counter1),
      m_counter2(counter2)
    {}

    /*! Operate */
    Arr_parameter_space operator()(const Point_2& p) const
    { ++m_counter1; return m_object(p); }


    /*! Operate */
    Arr_parameter_space operator()(const X_monotone_curve_2& xc) const
    { ++m_counter2; return m_object(xc); }
  };

  /*! A functor that compares the x-coordinate of two given points
   * that lie on horizontal boundaries.
   */
  class Compare_x_on_boundary_2 {
  private:
    typename Base::Compare_x_on_boundary_2 m_object;
    size_t& m_counter1;
    size_t& m_counter2;
    size_t& m_counter3;

  public:
    /*! Construct */
  Compare_x_on_boundary_2(const Base* base,  size_t& counter1,
                          size_t& counter2, size_t& counter3 ) :
      m_object(base->compare_x_on_boundary_2_object()),
      m_counter1(counter1),
      m_counter2(counter2),
      m_counter3(counter3)
    {}

    /*! Operate */
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const
    { ++m_counter1; return m_object(p1, p2); }

    /*! Operate */
    Comparison_result operator()(const Point_2& pt,
                                 const X_monotone_curve_2& xcv,
                                 Arr_curve_end ce) const
    { ++m_counter2; return m_object(pt, xcv, ce); }

    /*! Operate */
    Comparison_result operator()(const X_monotone_curve_2& xcv1,
                                 Arr_curve_end ce1,
                                 const X_monotone_curve_2& xcv2,
                                 Arr_curve_end ce2) const
    { ++m_counter3; return m_object(xcv1, ce1, xcv2, ce2); }
  };

  /*! A functor that compares the x-coordinates of curve ends near the
   * boundary of the parameter space.
   */
  class Compare_x_near_boundary_2 {
  private:
    typename Base::Compare_x_near_boundary_2 m_object;
    size_t& m_counter;

  public:
    /*! Construct */
    Compare_x_near_boundary_2(const Base* base, size_t& counter) :
      m_object(base->compare_x_near_boundary_2_object()),
      m_counter(counter)
    {}

    /*! Operate */
    Comparison_result operator()(const X_monotone_curve_2& xc1,
                                 const X_monotone_curve_2& xc2,
                                 Arr_curve_end ce) const
    { ++m_counter; return m_object(xc1, xc2, ce); }
  };

  //@}



  /// \name Obtain the appropriate functor
  //@{

  Compare_x_2 compare_x_2_object() const
  { return Compare_x_2(this, m_counters[COMPARE_X_OP]); }

  Compare_xy_2 compare_xy_2_object() const
  { return Compare_xy_2(this, m_counters[COMPARE_XY_OP]); }

  Construct_min_vertex_2 construct_min_vertex_2_object() const
  { return Construct_min_vertex_2(this, m_counters[CONSTRUCT_MIN_VERTEX_OP]); }

  Construct_max_vertex_2 construct_max_vertex_2_object() const
  { return Construct_max_vertex_2(this, m_counters[CONSTRUCT_MAX_VERTEX_OP]); }

  Is_vertical_2 is_vertical_2_object() const
  { return Is_vertical_2(this, m_counters[IS_VERTICAL_OP]); }

  Compare_y_at_x_2 compare_y_at_x_2_object() const
  { return Compare_y_at_x_2(this, m_counters[COMPARE_Y_AT_X_OP]); }

  Equal_2 equal_2_object() const
  {
    return Equal_2(this, m_counters[EQUAL_POINTS_OP],
                   m_counters[EQUAL_CURVES_OP]);
  }

  Compare_y_at_x_left_2 compare_y_at_x_left_2_object() const
  { return Compare_y_at_x_left_2(this, m_counters[COMPARE_Y_AT_X_LEFT_OP]); }

  Compare_y_at_x_right_2 compare_y_at_x_right_2_object() const
  { return Compare_y_at_x_right_2(this, m_counters[COMPARE_Y_AT_X_RIGHT_OP]); }

  Make_x_monotone_2 make_x_monotone_2_object() const
  { return Make_x_monotone_2(this, m_counters[MAKE_X_MONOTONE_OP]); }

  Split_2 split_2_object() const
  { return Split_2(this, m_counters[SPLIT_OP]); }

  Intersect_2 intersect_2_object() const
  { return Intersect_2(this, m_counters[INTERSECT_OP]); }

  Are_mergeable_2 are_mergeable_2_object() const
  { return Are_mergeable_2(this, m_counters[ARE_MERGEABLE_OP]); }

  Merge_2 merge_2_object() const
  { return Merge_2(this, m_counters[MERGE_OP]); }

  Construct_opposite_2 construct_opposite_2_object() const
  { return Construct_opposite_2(this, m_counters[CONSTRUCT_OPPOSITE_OP]); }

  Compare_endpoints_xy_2 compare_endpoints_xy_2_object() const
  { return Compare_endpoints_xy_2(this, m_counters[COMPARE_ENDPOINTS_XY_OP]); }

  // left-right
  Parameter_space_in_x_2 parameter_space_in_x_2_object() const
  { return Parameter_space_in_x_2(
        this,
        m_counters[PARAMETER_SPACE_IN_X_CURVE_END_OP],
        m_counters[PARAMETER_SPACE_IN_X_POINT_OP],
        m_counters[PARAMETER_SPACE_IN_X_CURVE_OP]
    );
  }

  Is_on_x_identification_2 is_on_x_identification_2_object() const
  {
    return Is_on_x_identification_2(this,
                                    m_counters[IS_ON_X_IDENTIFICATION_POINT_OP],
                                    m_counters[IS_ON_X_IDENTIFICATION_CURVE_OP]);
  }

  Compare_y_on_boundary_2 compare_on_boundary_2_object() const
  { return Compare_y_on_boundary_2(this, m_counters[COMPARE_Y_ON_BOUNDARY_OP]); }

  Compare_y_near_boundary_2 compare_near_boundary_2_object() const
  {
    return Compare_y_near_boundary_2(this,
                                     m_counters[COMPARE_Y_NEAR_BOUNDARY_OP]);
  }

  // bottom-top
  Parameter_space_in_y_2 parameter_space_in_y_2_object() const
  { return Parameter_space_in_y_2(
        this,
        m_counters[PARAMETER_SPACE_IN_Y_CURVE_END_OP],
        m_counters[PARAMETER_SPACE_IN_Y_POINT_OP],
        m_counters[PARAMETER_SPACE_IN_Y_CURVE_OP]
    );
  }

  Is_on_y_identification_2 is_on_y_identification_2_object() const
  { return Is_on_y_identification_2(
        this,
        m_counters[IS_ON_Y_IDENTIFICATION_POINT_OP],
        m_counters[IS_ON_Y_IDENTIFICATION_CURVE_OP]
    );
  }

  Compare_x_on_boundary_2 compare_x_on_boundary_2_object() const
  {
    return
      Compare_x_on_boundary_2(this,
                              m_counters[COMPARE_X_ON_BOUNDARY_POINTS_OP],
                              m_counters[COMPARE_X_ON_BOUNDARY_POINT_CURVE_END_OP],
                              m_counters[COMPARE_X_ON_BOUNDARY_CURVE_ENDS_OP]);
  }

  Compare_x_near_boundary_2 compare_x_near_boundary_2_object() const
  {
    return Compare_x_near_boundary_2(this,
                                     m_counters[COMPARE_X_NEAR_BOUNDARY_OP]);
  }

  //@}

  /*! Increment the construction counter
   * \param doit indicates whethet to actually inceremnt the counter or not
   * \return the counter at the end of the operation
   */
  static size_t increment(bool doit = true)
  {
#ifdef CGAL_NO_ATOMIC
    static size_t counter;
#else
    static std::atomic<size_t> counter;
#endif
    if (doit) ++counter;
    return counter;
  }

  /*! Clean all operation counters */
  void clear_counters()
  { m_counters = {}; }

private:
  /*! The operation counters */
  mutable std::array<size_t, NUMBER_OF_OPERATIONS> m_counters;
};

template <typename Out_stream, class Base_traits>
inline
Out_stream& operator<<(Out_stream& os,
                       const Arr_counting_traits_2<Base_traits>& traits)
{
  typedef Arr_counting_traits_2<Base_traits>            Traits;
  size_t sum = 0;
  size_t i;
  for (i = 0; i < Traits::NUMBER_OF_OPERATIONS; ++i)
    sum += traits.count(static_cast<typename Traits::Operation_id>(i));
  os << "# of COMPARE_X operation = "
     << traits.count_compare_x() << std::endl
     << "# of COMPARE_XY operation = "
     << traits.count_compare_xy() << std::endl
     << "# of CONSTRUCT_MIN_VERTEX operation = "
     << traits.count_construct_min_vertex() << std::endl
     << "# of CONSTRUCT_MAX_VERTEX operation = "
     << traits.count_construct_max_vertex() << std::endl
     << "# of IS_VERTICAL operation = "
     << traits.count_is_vertical() << std::endl
     << "# of COMPARE_Y_AT_X operation = "
     << traits.count_compare_y_at_x() << std::endl
     << "# of EQUAL_POINTS operation = "
     << traits.count_equal_points() << std::endl
     << "# of EQUAL_CURVES operation = "
     << traits.count_equal_curves() << std::endl
     << "# of COMPARE_Y_AT_X_LEFT operation = "
     << traits.count_compare_y_at_x_left() << std::endl
     << "# of COMPARE_Y_AT_X_RIGHT operation = "
     << traits.count_compare_y_at_x_right() << std::endl
     << "# of MAKE_X_MONOTONE operation = "
     << traits.count_make_x_monotone() << std::endl
     << "# of SPLIT operation = "
     << traits.count_split() << std::endl
     << "# of INTERSECT operation = "
     << traits.count_intersect() << std::endl
     << "# of ARE_MERGEABLE operation = "
     << traits.count_are_mergeable() << std::endl
     << "# of MERGE operation = "
     << traits.count_merge() << std::endl
     << "# of CONSTRUCT_OPPOSITE operation = "
     << traits.count_construct_opposite() << std::endl
     << "# of COMPARE_ENDPOINTS_XY operation = "
     << traits.count_compare_endpoints_xy() << std::endl
    // left-right
     << "# of PARAMETER_SPACE_IN_X curve-end operation = "
     << traits.count_parameter_space_in_x_curve_end() << std::endl
     << "# of PARAMETER_SPACE_IN_X point operation = "
     << traits.count_parameter_space_in_x_point() << std::endl
     << "# of PARAMETER_SPACE_IN_X curve operation = "
     << traits.count_parameter_space_in_x_curve() << std::endl
     << "# of IS_ON_X_IDENTIFICIATION point operation = "
     << traits.count_is_on_x_identification_point() << std::endl
     << "# of IS_ON_X_IDENTIFICATION curve operation = "
     << traits.count_is_on_x_identification_curve() << std::endl
     << "# of COMPARE_Y_ON_BOUNDARY operation = "
     << traits.count_compare_y_on_boundary() << std::endl
     << "# of COMPARE_Y_NEAR_BOUNDARY operation = "
     << traits.count_compare_y_near_boundary() << std::endl
    // bottom-top
     << "# of PARAMETER_SPACE_IN_Y curve-end operation = "
     << traits.count_parameter_space_in_y_curve_end() << std::endl
     << "# of PARAMETER_SPACE_IN_Y point operation = "
     << traits.count_parameter_space_in_y_point() << std::endl
     << "# of PARAMETER_SPACE_IN_Y curve operation = "
     << traits.count_parameter_space_in_y_curve() << std::endl
     << "# of IS_ON_Y_IDENTIFICIATION point operation = "
     << traits.count_is_on_y_identification_point() << std::endl
     << "# of IS_ON_Y_IDENTIFICATION curve operation = "
     << traits.count_is_on_y_identification_curve() << std::endl
     << "# of COMPARE_X_ON_BOUNDARY points operation = "
     << traits.count_compare_x_on_boundary_points() << std::endl
     << "# of COMPARE_X_ON_BOUNDARY point/curve-end operation = "
     << traits.count_compare_x_on_boundary_point_curve_end() << std::endl
     << "# of COMPARE_X_ON_BOUNDARY curve-ends operation = "
     << traits.count_compare_x_on_boundary_curve_ends() << std::endl
     << "# of COMPARE_X_NEAR_BOUNDARY operation = "
     << traits.count_compare_x_near_boundary() << std::endl

     << "total # = " << sum << std::endl
     << "# of traits constructed = " << Traits::increment(false)
     << std::endl;
  return os;
}

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
