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
 */

#include <iostream>
#include <string.h>
#include <atomic>
#include <array>

#include <CGAL/basic.h>
#include <CGAL/Arr_enums.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_has.h>

namespace CGAL {

/*! \class
 * A metadata traits-class decorator for the arrangement package. It counts the
 * number of invocations of traits-class functors. It is parameterized with
 * another traits class and inherits from it. For each traits method it
 * maintains a counter that counts the number of invocations into the method.
 *
 * It models all the concepts that the original traits models.
 */
template <typename BaseTraits>
class Arr_counting_traits_2 : public BaseTraits {
public:
  enum Operation_id {
    COMPARE_X_2_OP = 0,
    COMPARE_XY_2_OP,
    CONSTRUCT_MIN_VERTEX_2_OP,
    CONSTRUCT_MAX_VERTEX_2_OP,
    IS_VERTICAL_2_OP,
    COMPARE_Y_AT_X_2_OP,
    EQUAL_2_POINTS_OP,
    EQUAL_2_CURVES_OP,
    COMPARE_Y_AT_X_LEFT_2_OP,
    COMPARE_Y_AT_X_RIGHT_2_OP,
    MAKE_X_MONOTONE_2_OP,
    SPLIT_2_OP,
    INTERSECT_2_OP,
    ARE_MERGEABLE_2_OP,
    MERGE_2_OP,
    CONSTRUCT_2_OPPOSITE_2_OP,
    COMPARE_ENDPOINTS_XY_2_OP,
    APPROXIMATE_2_COORD_OP,
    APPROXIMATE_2_POINT_OP,
    APPROXIMATE_2_CURVE_OP,

    PARAMETER_SPACE_IN_X_2_CURVE_END_OP,
    PARAMETER_SPACE_IN_X_2_POINT_OP,
    IS_ON_X_IDENTIFICATION_POINT_2_OP,
    IS_ON_X_IDENTIFICATION_CURVE_2_OP,
    COMPARE_Y_ON_BOUNDARY_2_OP,
    COMPARE_Y_NEAR_BOUNDARY_2_OP,

    PARAMETER_SPACE_IN_Y_2_CURVE_END_OP,
    PARAMETER_SPACE_IN_Y_2_POINT_OP,
    IS_ON_Y_IDENTIFICATION_2_POINT_OP,
    IS_ON_Y_IDENTIFICATION_2_CURVE_OP,
    COMPARE_X_ON_BOUNDARY_2_POINTS_OP,
    COMPARE_X_ON_BOUNDARY_2_POINT_CURVE_END_OP,
    COMPARE_X_ON_BOUNDARY_2_CURVE_ENDS_OP,
    COMPARE_X_NEAR_BOUNDARY_2_OP,

    NUMBER_OF_OPERATIONS
  };

  using Base = BaseTraits;

  /*! constructs default */
  template <typename ... Args>
  Arr_counting_traits_2(Args ... args) : Base(std::forward<Args>(args)...) {
    clear_counters();
    increment();
  }

  /*! disables copy constructor. */
  Arr_counting_traits_2(const Arr_counting_traits_2&) = delete;

  /*! obtains the counter of the given operation */
  std::size_t count(Operation_id id) const { return m_counters[id]; }

  /*! prints the counter associated with an operation. */
  template <typename OutStream>
  OutStream& print(OutStream& os, Operation_id id) const {
    if (! m_exist[id]) return os;
    os << m_names[id] << ": " << m_counters[id] << std::endl;
    return os;
  }

  /// \name Types and functors inherited from `BaseTraits`
  //@{

  // Traits types:
  using Has_left_category = typename Base::Has_left_category;
  using Has_merge_category = typename Base::Has_merge_category;
  using Has_do_intersect_category = typename Base::Has_do_intersect_category;

  using Left_side_category =
    typename internal::Arr_complete_left_side_category<Base>::Category;
  using Bottom_side_category =
    typename internal::Arr_complete_bottom_side_category<Base>::Category;
  using Top_side_category =
    typename internal::Arr_complete_top_side_category<Base>::Category;
  using Right_side_category =
    typename internal::Arr_complete_right_side_category<Base>::Category;

  using Point_2 = typename Base::Point_2;
  using X_monotone_curve_2 = typename Base::X_monotone_curve_2;
  using Curve_2 = typename Base::Curve_2;

  /*! A functor that compares the \f$x\f$-coordinates of two points */
  class Compare_x_2 {
  private:
    typename Base::Compare_x_2 m_object;
    std::size_t& m_counter;

  public:
    /*! constructs */
    Compare_x_2(const Base& base, std::size_t& counter) :
      m_object(base.compare_x_2_object()), m_counter(counter) {}

    /*! operates */
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const
    { ++m_counter; return m_object(p1, p2); }
  };

  /*! A functor that compares two points lexigoraphically: by \f$x\f$, then by
   * \f$y\f$. */
  class Compare_xy_2 {
  private:
    typename Base::Compare_xy_2 m_object;
    std::size_t& m_counter;

  public:
    /*! constructs */
    Compare_xy_2(const Base& base, std::size_t& counter) :
      m_object(base.compare_xy_2_object()), m_counter(counter) {}

    /*! operates */
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const
    { ++m_counter; return m_object(p1, p2); }
  };

  /*! A functor that obtains the left endpoint of an \f$x\f$-monotone curve. */
  class Construct_min_vertex_2 {
  private:
    typename Base::Construct_min_vertex_2 m_object;
    std::size_t& m_counter;

  public:
    /*! constructs */
    Construct_min_vertex_2(const Base& base, std::size_t& counter) :
      m_object(base.construct_min_vertex_2_object()), m_counter(counter) {}

    /*! operates */
    using Subcurve_ctr_minv = typename Base::Construct_min_vertex_2;
    decltype(std::declval<Subcurve_ctr_minv>().
             operator()(std::declval<X_monotone_curve_2>()))
    operator()(const X_monotone_curve_2& xcv) const
    { ++m_counter; return m_object(xcv); }
  };

  /*! A functor that obtains the right endpoint of an \f$x\f$-monotone curve. */
  class Construct_max_vertex_2 {
  private:
    typename Base::Construct_max_vertex_2 m_object;
    std::size_t& m_counter;

  public:
    /*! constructs */
    Construct_max_vertex_2(const Base& base, std::size_t& counter) :
      m_object(base.construct_max_vertex_2_object()), m_counter(counter) {}

    /*! operates */
    using Subcurve_ctr_maxv = typename Base::Construct_max_vertex_2;
    decltype(std::declval<Subcurve_ctr_maxv>().
             operator()(std::declval<X_monotone_curve_2>()))
    operator()(const X_monotone_curve_2& xcv) const
    { ++m_counter; return m_object(xcv); }
  };

  /*! A functor that checks whether a given \f$x\f$-monotone curve is vertical.
   */
  class Is_vertical_2 {
  private:
    typename Base::Is_vertical_2 m_object;
    std::size_t& m_counter;

  public:
    /*! constructs */
    Is_vertical_2(const Base& base, std::size_t& counter) :
      m_object(base.is_vertical_2_object()), m_counter(counter) {}

    /*! operates */
    bool operator()(const X_monotone_curve_2& xc) const
    { ++m_counter; return m_object(xc); }
  };

  /*! A functor that compares the \f$y\f$-coordinates of a point and an
   * \f$x\f$-monotone curve at the point \f$x\f$-coordinate.
   */
  class Compare_y_at_x_2 {
  private:
    typename Base::Compare_y_at_x_2 m_object;
    std::size_t& m_counter;

  public:
    /*! constructs */
    Compare_y_at_x_2(const Base& base, std::size_t& counter) :
      m_object(base.compare_y_at_x_2_object()), m_counter(counter) {}

    /*! operates */
    Comparison_result operator()(const Point_2& p,
                                 const X_monotone_curve_2& xc) const
    { ++m_counter; return m_object(p, xc); }
  };

  /*! A functor that checks whether two points and two \f$x\f$-monotone curves
   * are identical.
   */
  class Equal_2 {
  private:
    typename Base::Equal_2 m_object;
    std::size_t& m_counter1;
    std::size_t& m_counter2;

  public:
    /*! constructs */
    Equal_2(const Base& base, std::size_t& counter1, std::size_t& counter2) :
      m_object(base.equal_2_object()),
      m_counter1(counter1), m_counter2(counter2)
    {}

    /*! operates */
    bool operator()(const Point_2& p1, const Point_2& p2) const
    { ++m_counter1; return m_object(p1, p2); }

    /*! operates */
    bool operator()(const X_monotone_curve_2& xc1,
                    const X_monotone_curve_2& xc2) const
    { ++m_counter2; return m_object(xc1, xc2); }
  };

  /*! A functor that compares compares the \f$y\f$-coordinates of two
   * \f$x\f$-monotone curves immediately to the left of their intersection
   * point.
   */
  class Compare_y_at_x_left_2 {
  private:
    typename Base::Compare_y_at_x_left_2 m_object;
    std::size_t& m_counter;

  public:
    /*! constructs */
    Compare_y_at_x_left_2(const Base& base, std::size_t& counter) :
      m_object(base.compare_y_at_x_left_2_object()), m_counter(counter) {}

    /*! operates */
    Comparison_result operator()(const X_monotone_curve_2& xc1,
                                 const X_monotone_curve_2& xc2,
                                 const Point_2& p) const
    { ++m_counter; return m_object(xc1, xc2, p); }
  };

  /*! A functor that compares compares the \f$y\f$-coordinates of two
   * \f$x\f$-monotone curves immediately to the right of their intersection
   * point.
   */
  class Compare_y_at_x_right_2 {
  private:
    typename Base::Compare_y_at_x_right_2 m_object;
    std::size_t& m_counter;

  public:
    /*! constructs */
    Compare_y_at_x_right_2(const Base& base, std::size_t& counter) :
      m_object(base.compare_y_at_x_right_2_object()), m_counter(counter) {}

    /*! operates */
    Comparison_result operator()(const X_monotone_curve_2& xc1,
                                 const X_monotone_curve_2& xc2,
                                 const Point_2& p) const
    { ++m_counter; return m_object(xc1, xc2, p); }
  };

  /*! \class Make_x_monotone_2
   * A functor that subdivides a curve into \f$x\f$-monotone curves.
   */
  class Make_x_monotone_2 {
  private:
    typename Base::Make_x_monotone_2 m_object;
    std::size_t& m_counter;

  public:
    /*! constructs */
    Make_x_monotone_2(const Base& base, std::size_t& counter) :
      m_object(base.make_x_monotone_2_object()), m_counter(counter) {}

    /*! subdivides a given curve into \f$x\f$-monotone subcurves and insert them
     * into a given output iterator.
     * \param cv the curve.
     * \param oi the output iterator for the result. Its value type is a variant
     *           that wraps `Point_2` or an `X_monotone_curve_2` objects.
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
    std::size_t& m_counter;

  public:
    /*! constructs */
    Split_2(const Base& base, std::size_t& counter) :
      m_object(base.split_2_object()), m_counter(counter) {}

    /*! operates */
    void operator()(const X_monotone_curve_2& xc, const Point_2& p,
                    X_monotone_curve_2& xc1, X_monotone_curve_2& xc2) const
    { ++m_counter; m_object(xc, p, xc1, xc2); }
  };

  /*! A functor that computes intersections between \f$x\f$-monotone curves. */
  class Intersect_2 {
  private:
    typename Base::Intersect_2 m_object;
    std::size_t& m_counter;

  public:
    /*! constructs */
    Intersect_2(const Base& base, std::size_t& counter) :
      m_object(base.intersect_2_object()), m_counter(counter) {}

    /*! operates */
    template <typename OutputIterator>
    OutputIterator operator()(const X_monotone_curve_2& xc1,
                              const X_monotone_curve_2& xc2,
                              OutputIterator oi) const
    { ++m_counter; return m_object(xc1, xc2, oi); }
  };

  /*! A functor that tests whether two \f$x\f$-monotone curves can be merged. */
  class Are_mergeable_2 {
  private:
    typename Base::Are_mergeable_2 m_object;
    std::size_t& m_counter;

  public:
    /*! constructs */
    Are_mergeable_2(const Base& base, std::size_t& counter) :
      m_object(base.are_mergeable_2_object()), m_counter(counter) {}

    /*! operates */
    bool operator()(const X_monotone_curve_2& xc1,
                    const X_monotone_curve_2& xc2) const
    { ++m_counter; return m_object(xc1, xc2); }
  };

  /*! A functor that merges two \f$x\f$-monotone curves into one. */
  class Merge_2 {
  private:
    typename Base::Merge_2 m_object;
    std::size_t& m_counter;

  public:
    /*! constructs */
    Merge_2(const Base& base, std::size_t& counter) :
      m_object(base.merge_2_object()), m_counter(counter) {}

    /*! operates */
    void operator()(const X_monotone_curve_2& xc1,
                    const X_monotone_curve_2& xc2,
                    X_monotone_curve_2& xc) const
    { ++m_counter; m_object(xc1, xc2, xc); }
  };

  /*! A fnuctor that constructs an opposite \f$x\f$-monotone curve. */
  class Construct_opposite_2 {
  private:
    typename Base::Construct_opposite_2 m_object;
    std::size_t& m_counter;

  public:
    /*! constructs */
    Construct_opposite_2(const Base& base, std::size_t& counter) :
      m_object(base.construct_opposite_2_object()), m_counter(counter) {}

    /*! operates */
    X_monotone_curve_2 operator()(const X_monotone_curve_2& xc)
    { ++m_counter; return m_object(xc); }
  };

  /*! A functor that compares the two endpoints of an \f$x\f$-monotone curve
   * lexigoraphically.
   */
  class Compare_endpoints_xy_2 {
  private:
    typename Base::Compare_endpoints_xy_2 m_object;
    std::size_t& m_counter;

  public:
    /*! constructs */
    Compare_endpoints_xy_2(const Base& base, std::size_t& counter) :
      m_object(base.compare_endpoints_xy_2_object()), m_counter(counter) {}

    /*! operates */
    Comparison_result operator()(const X_monotone_curve_2& xc)
    { ++m_counter; return m_object(xc); }
  };

  /*! A functor that approximates coordinates, points, and \f$x\f$-monotone
   * curves.
   */
  class Approximate_2 {
  private:
    using Approximate_number_type = typename Base::Approximate_number_type;
    using Approximate_point_2 = typename Base::Approximate_point_2;

    typename Base::Approximate_2 m_object;
    std::size_t& m_counter1;
    std::size_t& m_counter2;
    std::size_t& m_counter3;

  public:
    /*! constructs */
    Approximate_2(const Base& base, std::size_t& counter1,
                  std::size_t& counter2, std::size_t& counter3) :
      m_object(base.approximate_2_object()),
      m_counter1(counter1),
      m_counter2(counter2),
      m_counter3(counter3)
    {}

    /*! obtains an approximation of a point coordinate.
     * \param p the exact point.
     * \param i the coordinate index (either 0 or 1).
     * \pre `i` is either 0 or 1.
     * \return An approximation of `p`'s \f$x\f$-coordinate (if `i` == 0), or an
     *         approximation of `p`'s \f$y\f$-coordinate (if `i` == 1).
     */
    Approximate_number_type operator()(const Point_2& p, std::size_t i) const
    { ++m_counter1; return m_object(p, i); }

    /*! obtains an approximation of a point.
     */
    Approximate_point_2 operator()(const Point_2& p) const
    { ++m_counter2; return m_object(p); }

    /*! obtains an approximation of an \f$x\f$-monotone curve.
     */
    template <typename OutputIterator>
    OutputIterator operator()(const X_monotone_curve_2& xcv, double error,
                              OutputIterator oi, bool l2r = true) const
    { ++m_counter3; return m_object(xcv, error, oi, l2r); }
  };

  // left-right

  /*! A functor that determines whether an endpoint of an \f$x\f$-monotone curve
   * lies on a boundary of the parameter space along the \f$x\f$-axis.
   */
  class Parameter_space_in_x_2 {
  private:
    typename Base::Parameter_space_in_x_2 m_object;
    std::size_t& m_counter1;
    std::size_t& m_counter2;

  public:
    /*! constructs */
    Parameter_space_in_x_2(const Base& base,
                           std::size_t& counter1, std::size_t& counter2) :
      m_object(base.parameter_space_in_x_2_object()),
      m_counter1(counter1),
      m_counter2(counter2)
    {}

    /*! operates */
    Arr_parameter_space operator()(const X_monotone_curve_2& xc,
                                   Arr_curve_end ce) const
    { ++m_counter1; return m_object(xc, ce); }


    /*! operates */
    Arr_parameter_space operator()(const Point_2& p) const
    { ++m_counter2; return m_object(p); }
  };

  /*! A functor that determines whether a point or a curve lies on an
   * identification in x.
   */
  class Is_on_x_identification_2 {
  private:
    typename Base::Is_on_x_identification_2 m_object;
    std::size_t& m_counter1;
    std::size_t& m_counter2;

  public:
    /*! constructs */
    Is_on_x_identification_2(const Base& base,
                             std::size_t& counter1, std::size_t& counter2) :
      m_object(base.is_on_x_identification_2_object()),
      m_counter1(counter1),
      m_counter2(counter2)
    {}

    /*! operates */
    bool operator()(const Point_2& p) const
    { ++m_counter1; return m_object(p); }

    /*! operates */
    bool operator()(const X_monotone_curve_2& xc) const
    { ++m_counter2; return m_object(xc); }
  };

  /*! A functor that compares the \f$y\f$-coordinate of two given points
   * that lie on vertical boundaries.
   */
  class Compare_y_on_boundary_2 {
  private:
    typename Base::Compare_y_on_boundary_2 m_object;
    std::size_t& m_counter;

  public:
    /*! constructs */
    Compare_y_on_boundary_2(const Base& base, std::size_t& counter) :
      m_object(base.compare_y_on_boundary_2_object()), m_counter(counter) {}

    /*! operates */
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const
    { ++m_counter; return m_object(p1, p2); }
  };

  /*! A functor that compares the \f$y\f$-coordinates of curve ends near the
   * boundary of the parameter space.
   */
  class Compare_y_near_boundary_2 {
  private:
    typename Base::Compare_y_near_boundary_2 m_object;
    std::size_t& m_counter;

  public:
    /*! constructs */
    Compare_y_near_boundary_2(const Base& base, std::size_t& counter) :
      m_object(base.compare_y_near_boundary_2_object()), m_counter(counter) {}

    /*! operates */
    Comparison_result operator()(const X_monotone_curve_2& xc1,
                                 const X_monotone_curve_2& xc2,
                                 Arr_curve_end ce) const
    { ++m_counter; return m_object(xc1, xc2, ce); }
  };

  // bottom-top

  /*! A functor that determines whether an endpoint of an \f$x\f$-monotone arc
   * lies on a boundary of the parameter space along the \f$y\f$-axis.
   */
  class Parameter_space_in_y_2 {
  private:
    typename Base::Parameter_space_in_y_2 m_object;
    std::size_t& m_counter1;
    std::size_t& m_counter2;

  public:
    /*! constructs */
    Parameter_space_in_y_2(const Base& base,
                           std::size_t& counter1, std::size_t& counter2) :
      m_object(base.parameter_space_in_y_2_object()),
      m_counter1(counter1),
      m_counter2(counter2)
    {}

    /*! operates */
    Arr_parameter_space operator()(const X_monotone_curve_2& xc,
                                   Arr_curve_end ce) const
    { ++m_counter1; return m_object(xc, ce); }

    /*! operates */
    Arr_parameter_space operator()(const Point_2& p) const
    { ++m_counter2; return m_object(p); }
  };

  /*! A functor that determines whether a point or a curve lies on an
   * identification in \f$x\f$.
   */
  class Is_on_y_identification_2 {
  private:
    typename Base::Is_on_y_identification_2 m_object;
    std::size_t& m_counter1;
    std::size_t& m_counter2;

  public:
    /*! constructs */
    Is_on_y_identification_2(const Base& base,
                             std::size_t& counter1, std::size_t& counter2) :
      m_object(base.is_on_y_identification_2_object()),
      m_counter1(counter1),
      m_counter2(counter2)
    {}

    /*! operates */
    bool operator()(const Point_2& p) const
    { ++m_counter1; return m_object(p); }


    /*! operates */
    bool operator()(const X_monotone_curve_2& xc) const
    { ++m_counter2; return m_object(xc); }
  };

  /*! A functor that compares the \f$x\f$-coordinate of two given points
   * that lie on horizontal boundaries.
   */
  class Compare_x_on_boundary_2 {
  private:
    typename Base::Compare_x_on_boundary_2 m_object;
    std::size_t& m_counter1;
    std::size_t& m_counter2;
    std::size_t& m_counter3;

  public:
    /*! constructs */
    Compare_x_on_boundary_2(const Base& base,  std::size_t& counter1,
                            std::size_t& counter2, std::size_t& counter3 ) :
      m_object(base.compare_x_on_boundary_2_object()),
      m_counter1(counter1),
      m_counter2(counter2),
      m_counter3(counter3)
    {}

    /*! operates */
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const
    { ++m_counter1; return m_object(p1, p2); }

    /*! operates */
    Comparison_result operator()(const Point_2& pt,
                                 const X_monotone_curve_2& xcv,
                                 Arr_curve_end ce) const
    { ++m_counter2; return m_object(pt, xcv, ce); }

    /*! operates */
    Comparison_result operator()(const X_monotone_curve_2& xcv1,
                                 Arr_curve_end ce1,
                                 const X_monotone_curve_2& xcv2,
                                 Arr_curve_end ce2) const
    { ++m_counter3; return m_object(xcv1, ce1, xcv2, ce2); }
  };

  /*! A functor that compares the \f$x\f$-coordinates of curve ends near the
   * boundary of the parameter space.
   */
  class Compare_x_near_boundary_2 {
  private:
    typename Base::Compare_x_near_boundary_2 m_object;
    std::size_t& m_counter;

  public:
    /*! constructs */
    Compare_x_near_boundary_2(const Base& base, std::size_t& counter) :
      m_object(base.compare_x_near_boundary_2_object()),
      m_counter(counter)
    {}

    /*! operates */
    Comparison_result operator()(const X_monotone_curve_2& xc1,
                                 const X_monotone_curve_2& xc2,
                                 Arr_curve_end ce) const
    { ++m_counter; return m_object(xc1, xc2, ce); }
  };

  //@}

  /// \name Obtain the appropriate functor
  //@{

  Compare_x_2 compare_x_2_object() const
  { return Compare_x_2(*this, m_counters[COMPARE_X_2_OP]); }

  Compare_xy_2 compare_xy_2_object() const
  { return Compare_xy_2(*this, m_counters[COMPARE_XY_2_OP]); }

  Construct_min_vertex_2 construct_min_vertex_2_object() const
  { return Construct_min_vertex_2(*this, m_counters[CONSTRUCT_MIN_VERTEX_2_OP]); }

  Construct_max_vertex_2 construct_max_vertex_2_object() const
  { return Construct_max_vertex_2(*this, m_counters[CONSTRUCT_MAX_VERTEX_2_OP]); }

  Is_vertical_2 is_vertical_2_object() const
  { return Is_vertical_2(*this, m_counters[IS_VERTICAL_2_OP]); }

  Compare_y_at_x_2 compare_y_at_x_2_object() const
  { return Compare_y_at_x_2(*this, m_counters[COMPARE_Y_AT_X_2_OP]); }

  Equal_2 equal_2_object() const {
    return Equal_2(*this, m_counters[EQUAL_2_POINTS_OP],
                   m_counters[EQUAL_2_CURVES_OP]);
  }

  Compare_y_at_x_left_2 compare_y_at_x_left_2_object() const
  { return Compare_y_at_x_left_2(*this, m_counters[COMPARE_Y_AT_X_LEFT_2_OP]); }

  Compare_y_at_x_right_2 compare_y_at_x_right_2_object() const
  { return Compare_y_at_x_right_2(*this, m_counters[COMPARE_Y_AT_X_RIGHT_2_OP]); }

  Make_x_monotone_2 make_x_monotone_2_object() const
  { return Make_x_monotone_2(*this, m_counters[MAKE_X_MONOTONE_2_OP]); }

  Split_2 split_2_object() const
  { return Split_2(*this, m_counters[SPLIT_2_OP]); }

  Intersect_2 intersect_2_object() const
  { return Intersect_2(*this, m_counters[INTERSECT_2_OP]); }

  Are_mergeable_2 are_mergeable_2_object() const
  { return Are_mergeable_2(*this, m_counters[ARE_MERGEABLE_2_OP]); }

  Merge_2 merge_2_object() const
  { return Merge_2(*this, m_counters[MERGE_2_OP]); }

  Construct_opposite_2 construct_opposite_2_object() const
  { return Construct_opposite_2(*this, m_counters[CONSTRUCT_2_OPPOSITE_2_OP]); }

  Compare_endpoints_xy_2 compare_endpoints_xy_2_object() const
  { return Compare_endpoints_xy_2(*this, m_counters[COMPARE_ENDPOINTS_XY_2_OP]); }

  Approximate_2 approximate_2_object() const {
    return Approximate_2(*this,
                         m_counters[APPROXIMATE_2_COORD_OP],
                         m_counters[APPROXIMATE_2_POINT_OP],
                         m_counters[APPROXIMATE_2_CURVE_OP]);
  }

  // left-right
  Parameter_space_in_x_2 parameter_space_in_x_2_object() const {
    return Parameter_space_in_x_2(*this,
                                  m_counters[PARAMETER_SPACE_IN_X_2_CURVE_END_OP],
                                  m_counters[PARAMETER_SPACE_IN_X_2_POINT_OP]);
  }

  Is_on_x_identification_2 is_on_x_identification_2_object() const {
    return Is_on_x_identification_2(*this,
                                    m_counters[IS_ON_X_IDENTIFICATION_POINT_2_OP],
                                    m_counters[IS_ON_X_IDENTIFICATION_CURVE_2_OP]);
  }

  Compare_y_on_boundary_2 compare_y_on_boundary_2_object() const {
    return Compare_y_on_boundary_2(*this, m_counters[COMPARE_Y_ON_BOUNDARY_2_OP]);
  }

  Compare_y_near_boundary_2 compare_y_near_boundary_2_object() const {
    return Compare_y_near_boundary_2(*this,
                                     m_counters[COMPARE_Y_NEAR_BOUNDARY_2_OP]);
  }

  // bottom-top
  Parameter_space_in_y_2 parameter_space_in_y_2_object() const {
    return Parameter_space_in_y_2(*this,
                                  m_counters[PARAMETER_SPACE_IN_Y_2_CURVE_END_OP],
                                  m_counters[PARAMETER_SPACE_IN_Y_2_POINT_OP]);
  }

  Is_on_y_identification_2 is_on_y_identification_2_object() const {
    return Is_on_y_identification_2(*this,
                                    m_counters[IS_ON_Y_IDENTIFICATION_2_POINT_OP],
                                    m_counters[IS_ON_Y_IDENTIFICATION_2_CURVE_OP]);
  }

  Compare_x_on_boundary_2 compare_x_on_boundary_2_object() const {
    return
      Compare_x_on_boundary_2(*this,
                              m_counters[COMPARE_X_ON_BOUNDARY_2_POINTS_OP],
                              m_counters[COMPARE_X_ON_BOUNDARY_2_POINT_CURVE_END_OP],
                              m_counters[COMPARE_X_ON_BOUNDARY_2_CURVE_ENDS_OP]);
  }

  Compare_x_near_boundary_2 compare_x_near_boundary_2_object() const {
    return Compare_x_near_boundary_2(*this,
                                     m_counters[COMPARE_X_NEAR_BOUNDARY_2_OP]);
  }

  //@}

  /*! increments the construction counter.
   * \param doit indicates whethet to actually inceremnt the counter or not.
   * \return the counter at the end of the operation.
   */
  static std::size_t increment(bool doit = true) {
#ifdef CGAL_NO_ATOMIC
    static std::size_t counter;
#else
    static std::atomic<std::size_t> counter;
#endif
    if (doit) ++counter;
    return counter;
  }

  /*! cleans all operation counters */
  void clear_counters() { m_counters = {}; }

private:
  //! The operation counters
  mutable std::array<std::size_t, NUMBER_OF_OPERATIONS> m_counters;
  const std::array<std::string, NUMBER_OF_OPERATIONS> m_names = {
    "COMPARE_X_2_OP",
    "COMPARE_XY_2_OP",
    "CONSTRUCT_MIN_VERTEX_2_OP",
    "CONSTRUCT_MAX_VERTEX_2_OP",
    "IS_VERTICAL_2_OP",
    "COMPARE_Y_AT_X_2_OP",
    "EQUAL_2_POINTS_OP",
    "EQUAL_2_CURVES_OP",
    "COMPARE_Y_AT_X_LEFT_2_OP",
    "COMPARE_Y_AT_X_RIGHT_2_OP",
    "MAKE_X_MONOTONE_2_OP",
    "SPLIT_2_OP",
    "INTERSECT_2_OP",
    "ARE_MERGEABLE_2_OP",
    "MERGE_2_OP",
    "CONSTRUCT_2_OPPOSITE_2_OP",
    "COMPARE_ENDPOINTS_XY_2_OP",
    "APPROXIMATE_2_COORD_OP",
    "APPROXIMATE_2_POINT_OP",
    "APPROXIMATE_2_CURVE_OP",

    "PARAMETER_SPACE_IN_X_2_CURVE_END_OP",
    "PARAMETER_SPACE_IN_X_2_POINT_OP",
    "IS_ON_X_IDENTIFICATION_POINT_2_OP",
    "IS_ON_X_IDENTIFICATION_CURVE_2_OP",
    "COMPARE_Y_ON_BOUNDARY_2_OP",
    "COMPARE_Y_NEAR_BOUNDARY_2_OP",

    "PARAMETER_SPACE_IN_Y_2_CURVE_END_OP",
    "PARAMETER_SPACE_IN_Y_2_POINT_OP",
    "IS_ON_Y_IDENTIFICATION_2_POINT_OP",
    "IS_ON_Y_IDENTIFICATION_2_CURVE_OP",
    "COMPARE_X_ON_BOUNDARY_2_POINTS_OP",
    "COMPARE_X_ON_BOUNDARY_2_POINT_CURVE_END_OP",
    "COMPARE_X_ON_BOUNDARY_2_CURVE_ENDS_OP",
    "COMPARE_X_NEAR_BOUNDARY_2_OP"
  };
  const std::array<bool, NUMBER_OF_OPERATIONS> m_exist = {
    has_compare_x_2<Base>::value,
    has_compare_xy_2<Base>::value,
    has_construct_min_vertex_2<Base>::value,
    has_construct_max_vertex_2<Base>::value,
    has_is_vertical_2<Base>::value,
    has_compare_y_at_x_2<Base>::value,
    has_equal_2<Base>::value, // points
    has_equal_2<Base>::value, // curves
    has_compare_y_at_x_left_2<Base>::value,
    has_compare_y_at_x_right_2<Base>::value,
    has_make_x_monotone_2<Base>::value,
    has_split_2<Base>::value,
    has_intersect_2<Base>::value,
    has_are_mergeable_2<Base>::value,
    has_merge_2<Base>::value,
    has_construct_opposite_2<Base>::value,
    has_compare_endpoints_xy_2<Base>::value,
    has_approximate_2<Base>::value, // coordinate
    has_approximate_2<Base>::value, // point
    has_approximate_2<Base>::value, // curve
    has_parameter_space_in_x_2<Base>::value, // curve end
    has_parameter_space_in_x_2<Base>::value, // point
    has_is_on_x_identification_2<Base>::value, // point
    has_is_on_x_identification_2<Base>::value, // curve
    has_compare_y_on_boundary_2<Base>::value,
    has_compare_y_near_boundary_2<Base>::value,
    has_parameter_space_in_y_2<Base>::value, // curve end
    has_parameter_space_in_y_2<Base>::value, // point
    has_is_on_y_identification_2<Base>::value, // point
    has_is_on_y_identification_2<Base>::value, // curve
    has_compare_x_on_boundary_2<Base>::value, // points
    has_compare_x_on_boundary_2<Base>::value, // point, curve end
    has_compare_x_on_boundary_2<Base>::value, // curve ends
    has_compare_x_near_boundary_2<Base>::value
  };
};

template <typename OutStream, class BaseTraits>
inline OutStream& operator<<(OutStream& os,
                             const Arr_counting_traits_2<BaseTraits>& traits) {
  using Traits = Arr_counting_traits_2<BaseTraits>;
  std::size_t sum = 0;
  for (auto i = 0; i < Traits::NUMBER_OF_OPERATIONS; ++i) {
    sum += traits.count(static_cast<typename Traits::Operation_id>(i));
    traits.print(os, static_cast<typename Traits::Operation_id>(i));
  }
  os << "total # = " << sum << std::endl
     << "# of traits constructed = " << Traits::increment(false)
     << std::endl;
  return os;
}

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
