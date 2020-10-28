// Copyright (c) 2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Efi Fogel    <efif@post.tau.ac.il>

#ifndef CGAL_ARR_TRACING_TRAITS_H
#define CGAL_ARR_TRACING_TRAITS_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>

/*! \file
 * A tracing traits-class for the arrangement package.
 * This is a meta-traits class. It is parameterized with another traits class
 * and inherits from it. For each traits method it prints out its input
 * parameters and its output result
 */

#include <iostream>
#include <list>

#include <boost/variant.hpp>

#include <CGAL/basic.h>
#include <CGAL/Arr_enums.h>
#include <CGAL/Arr_tags.h>

namespace CGAL {

/*! \class
 * A model of the ArrangementTraits_2 concept that counts the methods invoked.
 */
template <typename Base_traits>
class Arr_tracing_traits_2 : public Base_traits {
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

    PARAMETER_SPACE_IN_X_OP,
    IS_ON_X_IDENTIFICATION_OP,
    COMPARE_Y_ON_BOUNDARY_OP,
    COMPARE_Y_NEAR_BOUNDARY_OP,

    PARAMETER_SPACE_IN_Y_OP,
    IS_ON_Y_IDENTIFICATION_OP,
    COMPARE_X_AT_LIMIT_OP,
    COMPARE_X_NEAR_LIMIT_OP,
    COMPARE_X_ON_BOUNDARY_OP,
    COMPARE_X_NEAR_BOUNDARY_OP,

    NUMBER_OF_OPERATIONS
  };

private:
  typedef Base_traits                           Base;
  typedef Arr_tracing_traits_2<Base>            Self;

  /*! A set of bits that indicate whether operations should be traced */
  unsigned int m_flags;

  bool compare_x_op() const
  { return (0 != (m_flags & (0x1 << COMPARE_X_OP))); }

  bool compare_xy_op() const
  { return (0 != (m_flags & (0x1 << COMPARE_XY_OP))); }

  bool construct_min_vertex_op() const
  { return (0 != (m_flags & (0x1 << CONSTRUCT_MIN_VERTEX_OP))); }

  bool construct_max_vertex_op() const
  { return (0 != (m_flags & (0x1 << CONSTRUCT_MAX_VERTEX_OP))); }

  bool is_vertical_op() const
  { return (0 != (m_flags & (0x1 << IS_VERTICAL_OP))); }

  bool compare_y_at_x_op() const
  { return (0 != (m_flags & (0x1 << COMPARE_Y_AT_X_OP))); }

  bool equal_points_op() const
  { return (0 != (m_flags & (0x1 << EQUAL_POINTS_OP))); }

  bool equal_curves_op() const
  { return (0 != (m_flags & (0x1 << EQUAL_CURVES_OP))); }

  bool compare_y_at_x_left_op() const
  { return (0 != (m_flags & (0x1 << COMPARE_Y_AT_X_LEFT_OP))); }

  bool compare_y_at_x_right_op() const
  { return (0 != (m_flags & (0x1 << COMPARE_Y_AT_X_RIGHT_OP))); }

  bool make_x_monotone_op() const
  { return (0 != (m_flags & (0x1 << MAKE_X_MONOTONE_OP))); }

  bool split_op() const
  { return (0 != (m_flags & (0x1 << SPLIT_OP))); }

  bool intersect_op() const
  { return (0 != (m_flags & (0x1 << INTERSECT_OP))); }

  bool are_mergeable_op() const
  { return (0 != (m_flags & (0x1 << ARE_MERGEABLE_OP))); }

  bool merge_op() const
  { return (0 != (m_flags & (0x1 << MERGE_OP))); }

  bool construct_opposite_op() const
  { return (0 != (m_flags & (0x1 << CONSTRUCT_OPPOSITE_OP))); }

  bool compare_endpoints_xy_op() const
  { return (0 != (m_flags & (0x1 << COMPARE_ENDPOINTS_XY_OP))); }

  // left-right

  bool parameter_space_in_x_op() const
  { return (0 != (m_flags & (0x1 << PARAMETER_SPACE_IN_X_OP))); }

  bool is_on_x_identification_op() const
  { return m_flags & (0x1 << IS_ON_X_IDENTIFICATION_OP); }

  bool compare_y_on_boundary_op() const
  { return (0 != (m_flags & (0x1 << COMPARE_Y_ON_BOUNDARY_OP))); }

  bool compare_y_near_boundary_op() const
  { return m_flags & (0x1 << COMPARE_Y_NEAR_BOUNDARY_OP); }

  // bottom-top

  bool parameter_space_in_y_op() const
  { return (0 != (m_flags & (0x1 << PARAMETER_SPACE_IN_Y_OP))); }

  bool is_on_y_identification_op() const
  { return m_flags & (0x1 << IS_ON_Y_IDENTIFICATION_OP); }

  bool compare_x_at_limit_op() const
  { return m_flags & (0x1 << COMPARE_X_AT_LIMIT_OP); }

  bool compare_x_near_limit_op() const
  { return m_flags & (0x1 << COMPARE_X_NEAR_LIMIT_OP); }

  bool compare_x_on_boundary_op() const
  { return (0 != (m_flags & (0x1 << COMPARE_X_ON_BOUNDARY_OP))); }

  bool compare_x_near_boundary_op() const
  { return m_flags & (0x1 << COMPARE_X_NEAR_BOUNDARY_OP); }

public:
  /*! Default constructor */
  Arr_tracing_traits_2() :
    Base()
  {
    enable_all_traces();
  }

  /*! Enable the trace of a traits operation
   * \param id the operation identifier
   */
  void enable_trace(Operation_id id) { m_flags |= 0x1 << id; }

  /*! Enable the trace of all traits operations
   * \param id the operation identifier
   */
  void enable_all_traces() { m_flags = 0xffffffff; }

  /*! Disable the trace of a traits operation
   * \param id the operation identifier
   */
  void disable_trace(Operation_id id) { m_flags &= ~(0x1 << id); }

  /*! Disable the trace of all traits operations
   * \param id the operation identifier
   */
  void disable_all_traces() { m_flags = 0x0; }

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

  typedef typename Base::Point_2                Point_2;
  typedef typename Base::X_monotone_curve_2     X_monotone_curve_2;
  typedef typename Base::Curve_2                Curve_2;
  typedef typename Base::Multiplicity           Multiplicity;

  /*! A functor that compares the x-coordinates of two points */
  class Compare_x_2 {
  private:
    typename Base::Compare_x_2 m_object;
    bool m_enabled;

  public:
    /*! Construct */
    Compare_x_2(const Base * base, bool enabled = true) :
      m_object(base->compare_x_2_object()), m_enabled(enabled) {}

    /*! Operate
     * \param p1 first point
     * \param p2 second point
     * \return the comparison result
     */
    Comparison_result operator()(const Point_2 & p1, const Point_2 & p2) const
    {
      if (!m_enabled) return m_object(p1, p2);
      std::cout << "compare_x" << std::endl
                << "  p1: " << p1 << std::endl
                << "  p2: " << p2 << std::endl;
      Comparison_result cr = m_object(p1, p2);
      std::cout << "  result: " << cr << std::endl;
      return cr;
    }
  };

  /*! A functor that compares two points lexigoraphically: by x, then by y. */
  class Compare_xy_2 {
  private:
    typename Base::Compare_xy_2 m_object;
    bool m_enabled;

  public:
    /*! Construct */
    Compare_xy_2(const Base * base, bool enabled = true) :
      m_object(base->compare_xy_2_object()), m_enabled(enabled) {}

    /*! Operate
     * \param p1 the first point
     * \param p2 the second point
     * \return the comparison result
     */
    Comparison_result operator()(const Point_2 & p1, const Point_2 & p2) const
    {
      if (!m_enabled) return m_object(p1, p2);
      std::cout << "compare_xy" << std::endl
                << "  p1: " << p1 << std::endl
                << "  p2: " << p2 << std::endl;
      Comparison_result cr = m_object(p1, p2);
      std::cout << "  result: " << cr << std::endl;
      return cr;
    }
  };

  /*! A functor that obtains the left endpoint of an x-monotone curve. */
  class Construct_min_vertex_2 {
  private:
    typename Base::Construct_min_vertex_2 m_object;
    bool m_enabled;

  public:
    /*! Construct */
    Construct_min_vertex_2(const Base * base, bool enabled = true) :
      m_object(base->construct_min_vertex_2_object()), m_enabled(enabled) {}

    /*! Operate
     * \param xcv the curev the left endpoint of which is obtained
     * \return the left endpoint
     */
    const Point_2 operator()(const X_monotone_curve_2 & xcv) const
    {
      if (!m_enabled) return m_object(xcv);
      std::cout << "construct_min_vertex" << std::endl
                << "  xcv: " << xcv << std::endl;
      Point_2 p = m_object(xcv);
      std::cout << "  result: " << p << std::endl;
      return p;
    }
  };

  /*! A functor that obtains the right endpoint of an x-monotone curve. */
  class Construct_max_vertex_2 {
  private:
    typename Base::Construct_max_vertex_2 m_object;
    bool m_enabled;

  public:
    /*! Construct */
    Construct_max_vertex_2(const Base * base, bool enabled = true) :
      m_object(base->construct_max_vertex_2_object()), m_enabled(enabled) {}

    /*! Operate
     * \param xcv the curev the right endpoint of which is obtained
     * \return the right endpoint
     */
    const Point_2 operator()(const X_monotone_curve_2 & xcv) const
    {
      if (!m_enabled) return m_object(xcv);
      std::cout << "construct_max_vertex" << std::endl
                << "  xcv: " << xcv << std::endl;
      Point_2 p = m_object(xcv);
      std::cout << "  result: " << p << std::endl;
      return p;
    }
  };

  /*! A functor that checks whether a given x-monotone curve is vertical. */
  class Is_vertical_2 {
  private:
    typename Base::Is_vertical_2 m_object;
    bool m_enabled;

  public:
    /*! Construct */
    Is_vertical_2(const Base * base, bool enabled = true) :
      m_object(base->is_vertical_2_object()), m_enabled(enabled) {}

    /*! Operate
     * \param xcv the curve
     * \return a Boolean that indicates whether the curve is vertical or not
     */
    bool operator()(const X_monotone_curve_2 & xcv) const
    {
      if (!m_enabled) return m_object(xcv);
      std::cout << "is_vertical" << std::endl
                << "  xcv: " << xcv << std::endl;
      bool is_vertical = m_object(xcv);
      std::cout << "  result: " << is_vertical << std::endl;
      return is_vertical;
    }
  };

  /*! A functor that compares the y-coordinates of a point and an
   * x-monotone curve at the point x-coordinate.
   */
  class Compare_y_at_x_2 {
  private:
    typename Base::Compare_y_at_x_2 m_object;
    bool m_enabled;

  public:
    /*! Construct */
    Compare_y_at_x_2(const Base * base, bool enabled = true) :
      m_object(base->compare_y_at_x_2_object()), m_enabled(enabled) {}

    /*! Operate
     * \param p the point
     * \param xcv the curve
     * \return the comparison result
     */
    Comparison_result operator()(const Point_2 & p,
                                 const X_monotone_curve_2 & xcv) const
    {
      if (!m_enabled) return m_object(p, xcv);
      std::cout << "compare_y_at_x" << std::endl
                << "  p: " << p << std::endl
                << "  xcv: " << xcv << std::endl;
      Comparison_result cr = m_object(p, xcv);
      std::cout << "  result: " << cr << std::endl;
      return cr;
    }
  };

  /*! A functor that checks whether two points and two x-monotone curves are
   * identical.
   */
  class Equal_2 {
  private:
    typename Base::Equal_2 m_object;
    bool m_enabled_point;
    bool m_enabled_curve;

  public:
    /*! Construct */
    Equal_2(const Base * base,
            bool enabled_point = true, bool enabled_curve = true) :
      m_object(base->equal_2_object()),
      m_enabled_point(enabled_point),
      m_enabled_curve(enabled_curve)
    {}

    /*! Operate
     * \param xcv1 the first curve
     * \param xcv2 the second curve
     * \return true if the x-monotone curves are equal and false otherwise
     */
    bool operator()(const X_monotone_curve_2 & xcv1,
                    const X_monotone_curve_2 & xcv2) const
    {
      if (!m_enabled_curve) return m_object(xcv1, xcv2);
      std::cout << "equal 1" << std::endl
                << "  xcv1: " << xcv1 << std::endl
                << "  xcv1: " << xcv1 << std::endl;
      bool equal = m_object(xcv1, xcv2);
      std::cout << "  result: " << equal << std::endl;
      return equal;
    }

    /*! Operate
     * \param p1 the first point
     * \param p2 the second point
     * \return true if the points are equal and false otherwise
     */
    bool operator()(const Point_2 & p1, const Point_2 & p2) const
    {
      if (!m_enabled_point) return m_object(p1, p2);
      std::cout << "equal 2" << std::endl
                << "  p1: " << p1 << std::endl
                << "  p2: " << p2 << std::endl;
      bool equal = m_object(p1, p2);
      std::cout << "  result: " << equal << std::endl;
      return equal;
    }
  };

  /*! A functor that compares compares the y-coordinates of two x-monotone
   * curves immediately to the left of their intersection point.
   */
  class Compare_y_at_x_left_2 {
  private:
    typename Base::Compare_y_at_x_left_2 m_object;
    bool m_enabled;

  public:
    /*! Construct */
    Compare_y_at_x_left_2(const Base * base, bool enabled = true) :
      m_object(base->compare_y_at_x_left_2_object()), m_enabled(enabled) {}

    /*! Operate
     * \param xcv1 the first curve
     * \param xcv2 the second curve
     * \param p the reference point
     * \return the comparison result
     */
    Comparison_result operator()(const X_monotone_curve_2 & xcv1,
                                 const X_monotone_curve_2 & xcv2,
                                 const Point_2 & p) const
    {
      if (!m_enabled) return m_object(xcv1, xcv2, p);
      std::cout << "compare_y_at_x_left" << std::endl
                << "  p: " << p << std::endl
                << "  xcv1: " << xcv1 << std::endl
                << "  xcv2: " << xcv2 << std::endl;
      Comparison_result cr = m_object(xcv1, xcv2, p);
      std::cout << "  result:" << cr << std::endl;
      return cr;
    }
  };

  /*! A functor that compares compares the y-coordinates of two x-monotone
   * curves immediately to the right of their intersection point.
   */
  class Compare_y_at_x_right_2 {
  private:
    typename Base::Compare_y_at_x_right_2 m_object;
    bool m_enabled;

  public:
    /*! Construct */
    Compare_y_at_x_right_2(const Base * base, bool enabled = true) :
      m_object(base->compare_y_at_x_right_2_object()), m_enabled(enabled) {}

    /*! Operate
     * \param xcv1 the first curve
     * \param xcv2 the second curve
     * \param p the reference point
     * \return the comparison result
     */
    Comparison_result operator()(const X_monotone_curve_2 & xcv1,
                                 const X_monotone_curve_2 & xcv2,
                                 const Point_2 & p) const
    {
      if (!m_enabled) return m_object(xcv1, xcv2, p);
      std::cout << "compare_y_at_x_right" << std::endl
                << "  p: " << p << std::endl
                << "  xcv1: " << xcv1 << std::endl
                << "  xcv2: " << xcv2 << std::endl;
      Comparison_result cr = m_object(xcv1, xcv2, p);
      std::cout << "  result: " << cr << std::endl;
      return cr;
    }
  };

  //! \name Intersections & subdivisions
  //@{

  /*! \class Make_x_monotone_2
   * A functor for subdividing curves into x-monotone curves.
   */
  class Make_x_monotone_2 {
  private:
    typename Base::Make_x_monotone_2 m_object;
    bool m_enabled;

  public:
    /*! Construct */
    Make_x_monotone_2(const Base * base, bool enabled = true) :
      m_object(base->make_x_monotone_2_object()), m_enabled(enabled) {}

    /*! Subdivide a given curve into x-monotone subcurves and insert them into
     * a given output iterator.
     * \param cv the curve.
     * \param oi an output iterator for the result. Its value type is a variant
     *           that wraps Point_2 or X_monotone_curve_2 objects.
     * \return the output iterator.
     */
    template <typename OutputIterator>
    OutputIterator operator()(const Curve_2 & cv, OutputIterator oi) const
    {
      if (! m_enabled) return m_object(cv, oi);
      std::cout << "make_x_monotone" << std::endl
                << "  cv: " << cv << std::endl;

      typedef boost::variant<Point_2, X_monotone_curve_2>
        Make_x_monotone_result;

      std::list<Make_x_monotone_result> container;
      m_object(cv, std::back_inserter(container));
      if (container.empty()) return oi;

      size_t i = 0;
      for (auto it = container.begin(); it != container.end(); ++it) {
        if (const auto* xcv = boost::get<X_monotone_curve_2>(*it)) {
          std::cout << "  result[" << i++ << "]: xcv: " << *xcv << std::endl;
          continue;
        }

        if (const Point_2* p = boost::get<Point_2>(*it)) {
          std::cout << "  result[" << i++ << "]: p: " << *p << std::endl;
          continue;
        }

        CGAL_error();
      }

      for (auto it = container.begin(); it != container.end(); ++it) *oi++ = *it;
      container.clear();
      return oi;
    }
  };

  /*! A functor that splits an x-monotone curve at a point. */
  class Split_2 {
  private:
    typename Base::Split_2 m_object;
    bool m_enabled;

  public:
    /*! Construct */
    Split_2(const Base * base, bool enabled = true) :
      m_object(base->split_2_object()), m_enabled(enabled) {}

    /*! Operate
     * \param xcv
     * \param p
     * \param xcv1
     * \param xcv2
     */
    void operator()(const X_monotone_curve_2 & xcv, const Point_2 & p,
                    X_monotone_curve_2 & xcv1, X_monotone_curve_2 & xcv2) const
    {
      if (!m_enabled) {
        m_object(xcv, p, xcv1, xcv2);
        return;
      }
      std::cout << "split: " << std::endl
                << "  xcv: " << xcv << std::endl
                << "  p: " << p << std::endl;
      m_object(xcv, p, xcv1, xcv2);
      std::cout << "  result xcv1: " << xcv1 << std::endl
                << "         xcv2: " << xcv2 << std::endl;
    }
  };

  /*! A functor that computes intersections between two x-monotone curves. */
  class Intersect_2 {
  private:
    typename Base::Intersect_2 m_object;
    bool m_enabled;

  public:
    /*! Construct */
    Intersect_2(const Base* base, bool enabled = true) :
      m_object(base->intersect_2_object()), m_enabled(enabled) {}

    /*! Compute the intersections of the two given curves and insert them into
     * a given output iterator.
     * \param xcv1 the first curve
     * \param xcv2 the ssecond curve
     * \param oi the output iterator for the result. It value type is a variant
     *           that wraps an x-monotone overlapping curve or a pair that
     *           consists of the intersection point and its multiplicity
     * \return the past-the-end output iterator.
     */
    template <typename OutputIterator>
    OutputIterator operator()(const X_monotone_curve_2 & xcv1,
                              const X_monotone_curve_2 & xcv2,
                              OutputIterator oi) const
    {
      typedef std::pair<Point_2, Multiplicity>          Intersection_point;
      typedef boost::variant<Intersection_point, X_monotone_curve_2>
                                                        Intersection_result;

      if (! m_enabled) return m_object(xcv1, xcv2, oi);

      std::cout << "intersect" << std::endl
                << "  xcv1: " << xcv1 << std::endl
                << "  xcv2: " << xcv2 << std::endl;
      std::list<Intersection_result> container;
      m_object(xcv1, xcv2, std::back_inserter(container));
      if (container.empty()) return oi;

      unsigned int i = 0;
      for (const auto& item : container) {
        const X_monotone_curve_2* xcv = boost::get<X_monotone_curve_2>(&item);
        if (xcv != nullptr) {
          std::cout << "  result[" << i++ << "]: xcv: " << *xcv << std::endl;
          continue;
        }

        const Intersection_point* ip = boost::get<Intersection_point>(&item);
        if (ip != nullptr) {
          std::cout << "  result[" << i++ << "]: p: " << ip->first
                    << ", multiplicity: " << ip->second << std::endl;
          continue;
        }
      }

      for (auto it = container.begin(); it != container.end(); ++it) *oi++ = *it;
      container.clear();
      return oi;
    }
  };

  /*! A functor that tests whether two x-monotone curves can be merged. */
  class Are_mergeable_2 {
  private:
    typename Base::Are_mergeable_2 m_object;
    bool m_enabled;

  public:
    /*! Construct */
    Are_mergeable_2(const Base * base, bool enabled = true) :
      m_object(base->are_mergeable_2_object()), m_enabled(enabled) {}

    /*! Operate
     * \param xcv1 the first curve
     * \param xcv2 the second curve
     * \return true if the the two curve are mergeable and false otherwise.
     * Two curves are mergeable if they have the same underlying theoretical
     * curve
     */
    bool operator()(const X_monotone_curve_2 & xcv1,
                    const X_monotone_curve_2 & xcv2) const
    {
      if (!m_enabled) return m_object(xcv1, xcv2);
      std::cout << "are_mergeable" << std::endl
                << "  xcv1: " << xcv1 << std::endl
                << "  xcv2: " << xcv2 << std::endl;
      bool are_mergeable = m_object(xcv1, xcv2);
      std::cout << "  result: " << are_mergeable << std::endl;
      return are_mergeable;
    }
  };

  /*! A functor that merges two x-monotone curves into one. */
  class Merge_2 {
  private:
    typename Base::Merge_2 m_object;
    bool m_enabled;

  public:
    /*! Construct */
    Merge_2(const Base * base, bool enabled = true) :
      m_object(base->merge_2_object()), m_enabled(enabled) {}

    /*! Operate
     * \param xcv1 the first curve
     * \param xcv2 the second curve
     * \param xcv the merged curve
     */
    void operator()(const X_monotone_curve_2 & xcv1,
                    const X_monotone_curve_2 & xcv2,
                    X_monotone_curve_2 & xcv) const
    {
      std::cout << "merge" << std::endl
                << "  xcv1: " << xcv1 << std::endl
                << "  xcv2: " << xcv2 << std::endl;
      return m_object(xcv1, xcv2, xcv);
      std::cout << "  result: " << xcv << std::endl;
    }
  };

  /*! A fnuctor that constructs an opposite x-monotone curve. */
  class Construct_opposite_2 {
  private:
    typename Base::Construct_opposite_2 m_object;
    bool m_enabled;

  public:
    /*! Construct */
    Construct_opposite_2(const Base * base, bool enabled = true) :
      m_object(base->construct_opposite_2_object()), m_enabled(enabled) {}

    /*! Operate
     * \param xcv the curve
     * \return the opposite curve
     */
    X_monotone_curve_2 operator()(const X_monotone_curve_2 & xcv)
    {
      if (!m_enabled) return m_object(xcv);
      std::cout << "construct_opposite" << std::endl
                << "  xcv: " << xcv << std::endl;
      X_monotone_curve_2 xcv_out = m_object(xcv);
      std::cout << "  result: " << xcv_out << std::endl;
      return xcv;
    }
  };

  /*! A functor that compares the two endpoints of an x-monotone curve
   * lexigoraphically.
   */
  class Compare_endpoints_xy_2 {
  private:
    typename Base::Compare_endpoints_xy_2 m_object;
    bool m_enabled;

  public:
    /*! Construct */
    Compare_endpoints_xy_2(const Base * base, bool enabled = true) :
      m_object(base->compare_endpoints_xy_2_object()), m_enabled(enabled) {}

    /*! Operate
     * \param xcv the curve
     * \return the comparison result
     */
    Comparison_result operator()(const X_monotone_curve_2 & xcv)
    {
      if (!m_enabled) return m_object(xcv);
      std::cout << "compare_endpoints_xy" << std::endl
                << "  xcv: " << xcv << std::endl;
      Comparison_result cr = m_object(xcv);
      std::cout << "  result: " << cr << std::endl;
      return cr;
    }
  };

  // left-right

  /*! A functor that determines whether an endpoint of an x-monotone curve lies
   * on a boundary of the parameter space along the x axis.
   */
  class Parameter_space_in_x_2 {
  private:
    typename Base::Parameter_space_in_x_2 m_object;
    bool m_enabled;

  public:
    /*! Construct */
    Parameter_space_in_x_2(const Base * base, bool enabled = true) :
      m_object(base->parameter_space_in_x_2_object()), m_enabled(enabled)
    {}

    /*! Operate
     * \param xcv the curve the end of which is tested
     * \param ce the curve-end identifier
     * \return the boundary type
     */
    Arr_parameter_space operator()(const X_monotone_curve_2 & xcv,
                             Arr_curve_end ce) const
    {
      if (!m_enabled) return m_object(xcv, ce);
      std::cout << "parameter_space_in_x" << std::endl
                << "  ce: " << ce << ", xcv: " << xcv << std::endl;
      Arr_parameter_space bt = m_object(xcv, ce);
      std::cout << "  result: " << bt << std::endl;
      return bt;
    }
  };

  /*! A functor that determines whether a point or curve is on
   * x-identification.
   */
  class Is_on_x_identification_2 {
  private:
    typename Base::Is_on_x_identification_2 m_object;
    bool m_enabled;

  public:
    /*! Construct */
    Is_on_x_identification_2(const Base * base, bool enabled = true) :
      m_object(base->is_on_x_identification_2_object()), m_enabled(enabled) {}
    /*! Operate
     * \param p1 the point.
     */
    Comparison_result operator()(const Point_2 & p) const
    {
      if (!m_enabled) return m_object(p);
      std::cout << "is_on_x_identification" << std::endl
                << "  p: " << p << std::endl;
      Comparison_result cr = m_object(p);
      std::cout << "  result: " << cr << std::endl;
      return cr;
    }

    /*! Operate
     * \param xcv1 the curve
     */
    Comparison_result operator()(const X_monotone_curve_2 & xcv) const
    {
      if (!m_enabled) return m_object(xcv);
      std::cout << "is_on_x_identification" << std::endl
                << "  xcv: " << xcv << std::endl;
      Comparison_result cr = m_object(xcv);
      std::cout << "  result: " << cr << std::endl;
      return cr;
    }
  };

  /*! A functor that compares the y-coordinate of two given points
   * that lie on vertical boundaries.
   */
  class Compare_y_on_boundary_2 {
  private:
    typename Base::Compare_y_on_boundary_2 m_object;
    bool m_enabled;

  public:
    /*! Construct */
    Compare_y_on_boundary_2(const Base * base, bool enabled = true) :
      m_object(base->compare_y_on_boundary_2_object()),
      m_enabled(enabled)
    {}

    /*! Operate
     * \param p1 the first point.
     * \param p2 the second point.
     */
    Comparison_result operator()(const Point_2 & p1, const Point_2 & p2) const
    {
      if (!m_enabled) return m_object(p1, p2);
      std::cout << "compare_y_on_boundary" << std::endl
                << "  p1: " << p1 << std::endl
                << "  p2: " << p2 << std::endl;
      Comparison_result cr = m_object(p1, p2);
      std::cout << "  result: " << cr << std::endl;
      return cr;
    }
  };

  /*! A functor that compares the y-coordinates of curve ends near the
   * boundary of the parameter space.
   */
  class Compare_y_near_boundary_2 {
  private:
    typename Base::Compare_y_near_boundary_2 m_object;
    bool m_enabled;

  public:
    /*! Construct */
    Compare_y_near_boundary_2(const Base * base, bool enabled = true) :
      m_object(base->compare_y_near_boundary_2_object()), m_enabled(enabled) {}

    /*! Operate
     * \param xcv1 the first curve the end point of which is tested
     * \param xcv2 the second curve the end point of which is tested
     * \param ce the curve-end identifier
     * \return the comparison result
     */
    Comparison_result operator()(const X_monotone_curve_2 & xcv1,
                                 const X_monotone_curve_2 & xcv2,
                                 Arr_curve_end ce) const
    {
      if (!m_enabled) return m_object(xcv1, xcv2, ce);
      std::cout << "compare_y_near_boundary" << std::endl
                << "  ce: " << ce << std::endl
                << "  xcv1: " << xcv1 << std::endl
                << "  xcv2: " << xcv2 << std::endl;
      Comparison_result cr = m_object(xcv1, xcv2, ce);
      std::cout << "  result: " << cr << std::endl;
      return cr;
    }
  };

  // bottom-top

  /*! A functor that determines whether an endpoint of an x-monotone arc lies
   * on a boundary of the parameter space along the y axis.
   */
  class Parameter_space_in_y_2 {
  private:
    typename Base::Parameter_space_in_y_2 m_object;
    bool m_enabled;

  public:
    /*! Construct */
    Parameter_space_in_y_2(const Base * base, bool enabled = true) :
      m_object(base->parameter_space_in_y_2_object()), m_enabled(enabled) {}

    /*! Operate
     * \param xcv the curve the end of which is tested
     * \param ce the curve-end identifier
     * \return the boundary type
     */
    Arr_parameter_space operator()(const X_monotone_curve_2 & xcv,
                             Arr_curve_end ce) const
    {
      if (!m_enabled) return m_object(xcv, ce);
        std::cout << "parameter_space_in_y" << std::endl
                  << "  ce: " << ce << ", xcv: " << xcv << std::endl;
      Arr_parameter_space bt = m_object(xcv, ce);
      std::cout << "  result: " << bt << std::endl;
      return bt;
    }

    /*! Operate
     * \param p the point
     * \return the boundary type
     */
    Arr_parameter_space operator()(const Point_2 & p) const
    {
      if (!m_enabled) return m_object(p);
        std::cout << "parameter_space_in_y" << std::endl
                  << "  point: " << p << std::endl;
      Arr_parameter_space bt = m_object(p);
      std::cout << "  result: " << bt << std::endl;
      return bt;
    }
  };

  /*! A functor that determines whether a point or curve is on
   * y-identification.
   */
  class Is_on_y_identification_2 {
  private:
    typename Base::Is_on_y_identification_2 m_object;
    bool m_enabled;

  public:
    /*! Construct */
    Is_on_y_identification_2(const Base * base, bool enabled = true) :
      m_object(base->is_on_y_identification_2_object()), m_enabled(enabled) {}
    /*! Operate
     * \param p1 the point.
     */
    Comparison_result operator()(const Point_2 & p) const
    {
      if (!m_enabled) return m_object(p);
      std::cout << "is_on_y_identification" << std::endl
                << "  p: " << p << std::endl;
      Comparison_result cr = m_object(p);
      std::cout << "  result: " << cr << std::endl;
      return cr;
    }

    /*! Operate
     * \param xcv1 the curve
     */
    Comparison_result operator()(const X_monotone_curve_2 & xcv) const
    {
      if (!m_enabled) return m_object(xcv);
      std::cout << "is_on_y_identification" << std::endl
                << "  xcv: " << xcv << std::endl;
      Comparison_result cr = m_object(xcv);
      std::cout << "  result: " << cr << std::endl;
      return cr;
    }
  };

  /*! A functor that compares the x-limits of curve ends on the
   * boundary of the parameter space.
   */
  class Compare_x_at_limit_2 {
  private:
    typename Base::Compare_x_at_limit_2 m_object;
    bool m_enabled;

  public:
    /*! Construct */
    Compare_x_at_limit_2(const Base * base, bool enabled = true) :
      m_object(base->compare_x_at_limit_2_object()), m_enabled(enabled) {}

    /*! Operate
     * \param p the first point
     * \param xcv the curve the end of which is to be compared
     * \param ce the curve-end identifier
     * \return the comparison result
     */
    Comparison_result operator()(const Point_2 & p,
                                 const X_monotone_curve_2 & xcv,
                                 Arr_curve_end ce) const
    {
      if (!m_enabled) return m_object(p, xcv, ce);
      std::cout << "compare_x_at_limit 1" << std::endl
                << "  p: " << p << std::endl
                << "  ce: " << ce << ", xcv: " << xcv << std::endl;
      Comparison_result cr = m_object(p, xcv, ce);
      std::cout << "  result: " << std::endl;
      return cr;
    }

    /*! Operate
     * \param xcv1 the first curve the end of which is to be compared
     * \param ce1 the identifier of the end of the first curve
     * \param xcv2 the second curve the end of which is to be compared
     * \param ce2 the identifier of the end of the second curve
     * \return the comparison result
     */
    Comparison_result operator()(const X_monotone_curve_2 & xcv1,
                                 Arr_curve_end ce1,
                                 const X_monotone_curve_2 & xcv2,
                                 Arr_curve_end ce2) const
    {
      if (!m_enabled) return m_object(xcv1, ce1, xcv2, ce2);
      std::cout << "compare_x_at_limit 2" << std::endl
                << "  ce1: " << ce1 << ", xcv1: " << xcv1 << std::endl
                << "  ce2: " << ce2 << ", xcv2: " << xcv2 << std::endl;
      Comparison_result cr = m_object(xcv1, ce1, xcv2, ce2);
      std::cout << "  result: " << cr << std::endl;
      return cr;
    }
  };


  /*! A functor that compares the x-coordinates of curve ends near the
   * boundary of the parameter space.
   */
  class Compare_x_near_limit_2 {
  private:
    typename Base::Compare_x_near_limit_2 m_object;
    bool m_enabled;

  public:
    /*! Construct */
    Compare_x_near_limit_2(const Base * base, bool enabled = true) :
      m_object(base->compare_x_near_limit_2_object()), m_enabled(enabled) {}

    /*! Operate
     * \param xcv1 the first curve the end of which is to be compared
     * \param ce1 the identifier of the end of the first curve
     * \param xcv2 the second curve the end of which is to be compared
     * \param ce2 the identifier of the end of the second curve
     * \return the comparison result
     */
    Comparison_result operator()(const X_monotone_curve_2 & xcv1,
                                 const X_monotone_curve_2 & xcv2,
                                 Arr_curve_end ce) const
    {
      if (!m_enabled) return m_object(xcv1, xcv2, ce);
      std::cout << "compare_x_near_limit 2" << std::endl
                << "  xcv1: " << xcv1 << std::endl
                << "  xcv2: " << xcv2 << std::endl
                << "    ce: " << ce << std::endl;
      Comparison_result cr = m_object(xcv1, xcv2, ce);
      std::cout << "  result: " << cr << std::endl;
      return cr;
    }
  };

  /*! A functor that compares the x-coordinate of two given points
   * that lie on horizontal boundaries.
   */
  class Compare_x_on_boundary_2 {
  private:
    typename Base::Compare_x_on_boundary_2 m_object;
    bool m_enabled;

  public:
    /*! Construct */
    Compare_x_on_boundary_2(const Base * base, bool enabled = true) :
      m_object(base->compare_x_on_boundary_2_object()), m_enabled(enabled) {}
    /*! Operate
     * \param p1 the first point.
     * \param p2 the second point.
     */
    Comparison_result operator()(const Point_2 & p1, const Point_2 & p2) const
    {
      if (!m_enabled) return m_object(p1, p2);
      std::cout << "compare_x_on_boundary" << std::endl
                << "  p1: " << p1 << std::endl
                << "  p2: " << p2 << std::endl;
      Comparison_result cr = m_object(p1, p2);
      std::cout << "  result: " << cr << std::endl;
      return cr;
    }

    /*! Operate
     * \param pt the point.
     * \param xcv the curve.
     * \param ce the curve-end
     */
    Comparison_result operator()(const Point_2 & pt,
                                 const X_monotone_curve_2 & xcv, Arr_curve_end ce) const
    {
      if (!m_enabled) return m_object(pt, xcv, ce);
      std::cout << "compare_x_on_boundary" << std::endl
                << "  pt: " << pt << std::endl
                << " xcv: " << xcv << std::endl
                << "  ce: " << ce << std::endl;
      Comparison_result cr = m_object(pt, xcv, ce);
      std::cout << "  result: " << cr << std::endl;
      return cr;
    }

    /*! Operate
     * \param xcv1 the first curve.
     * \param ce1 the first curve-end
     * \param xcv2 the second curve.
     * \param ce2 the second curve-end
     */
    Comparison_result operator()(const X_monotone_curve_2 & xcv1, Arr_curve_end ce1,
                                 const X_monotone_curve_2 & xcv2, Arr_curve_end ce2) const
    {
      if (!m_enabled) return m_object(xcv2, ce1, xcv2, ce2);
      std::cout << "compare_x_on_boundary" << std::endl
                << "xcv1: " << xcv1 << std::endl
                << " ce1: " << ce1 << std::endl
                << "xcv2: " << xcv2 << std::endl
                << " ce2: " << ce2 << std::endl;
      Comparison_result cr = m_object(xcv1, ce1, xcv2, ce2);
      std::cout << "  result: " << cr << std::endl;
      return cr;
    }

  };

  /*! A functor that compares the x-coordinates of curve ends near the
   * boundary of the parameter space.
   */
  class Compare_x_near_boundary_2 {
  private:
    typename Base::Compare_x_near_boundary_2 m_object;
    bool m_enabled;

  public:
    /*! Construct */
    Compare_x_near_boundary_2(const Base * base, bool enabled = true) :
      m_object(base->compare_x_near_boundary_2_object()), m_enabled(enabled) {}

    /*! Operate
     * \param xcv1 the first curve the end of which is to be compared
     * \param ce1 the identifier of the end of the first curve
     * \param xcv2 the second curve the end of which is to be compared
     * \param ce2 the identifier of the end of the second curve
     * \return the comparison result
     */
    Comparison_result operator()(const X_monotone_curve_2 & xcv1,
                                 const X_monotone_curve_2 & xcv2,
                                 Arr_curve_end ce) const
    {
      if (!m_enabled) return m_object(xcv1, xcv2, ce);
      std::cout << "compare_x_near_boundary 2" << std::endl
                << "  xcv1: " << xcv1 << std::endl
                << "  xcv2: " << xcv2 << std::endl
                << "    ce: " << ce << std::endl;
      Comparison_result cr = m_object(xcv1, xcv2, ce);
      std::cout << "  result: " << cr << std::endl;
      return cr;
    }
  };

  //@}

  /// \name Obtain the appropriate functor
  //@{

  Compare_x_2 compare_x_2_object() const
  { return Compare_x_2(this, compare_x_op()); }

  Compare_xy_2 compare_xy_2_object() const
  { return Compare_xy_2(this, compare_xy_op()); }

  Construct_min_vertex_2 construct_min_vertex_2_object() const
  { return Construct_min_vertex_2(this, construct_min_vertex_op()); }

  Construct_max_vertex_2 construct_max_vertex_2_object() const
  { return Construct_max_vertex_2(this, construct_max_vertex_op()); }

  Is_vertical_2 is_vertical_2_object() const
  { return Is_vertical_2(this, is_vertical_op()); }

  Compare_y_at_x_2 compare_y_at_x_2_object() const
  { return Compare_y_at_x_2(this, compare_y_at_x_op()); }

  Equal_2 equal_2_object() const
  { return Equal_2(this, equal_points_op(), equal_curves_op()); }

  Compare_y_at_x_left_2 compare_y_at_x_left_2_object() const
  { return Compare_y_at_x_left_2(this, compare_y_at_x_left_op()); }

  Compare_y_at_x_right_2 compare_y_at_x_right_2_object() const
  { return Compare_y_at_x_right_2(this, compare_y_at_x_right_op()); }

  Make_x_monotone_2 make_x_monotone_2_object() const
  { return Make_x_monotone_2(this, make_x_monotone_op()); }

  Split_2 split_2_object() const
  { return Split_2(this, split_op()); }

  Intersect_2 intersect_2_object() const
  { return Intersect_2(this, intersect_op()); }

  Are_mergeable_2 are_mergeable_2_object() const
  { return Are_mergeable_2(this, are_mergeable_op()); }

  Merge_2 merge_2_object() const
  { return Merge_2(this, merge_op()); }

  Construct_opposite_2 construct_opposite_2_object() const
  { return Construct_opposite_2(this, construct_opposite_op()); }

  Compare_endpoints_xy_2 compare_endpoints_xy_2_object() const
  { return Compare_endpoints_xy_2(this, compare_endpoints_xy_op()); }

  // left-right

  Parameter_space_in_x_2 parameter_space_in_x_2_object() const
  { return Parameter_space_in_x_2(this, parameter_space_in_x_op()); }

  Is_on_x_identification_2 is_on_x_identification_2_object() const
  { return Is_on_x_identification_2(this, is_on_x_identification_op()); }

  Compare_y_on_boundary_2 compare_y_on_boundary_2_object() const
  { return Compare_y_on_boundary_2(this, compare_y_on_boundary_op()); }

  Compare_y_near_boundary_2 compare_y_near_boundary_2_object() const
  { return Compare_y_near_boundary_2(this, compare_y_near_boundary_op()); }

  // bottom-top

  Parameter_space_in_y_2 parameter_space_in_y_2_object() const
  { return Parameter_space_in_y_2(this, parameter_space_in_y_op()); }

  Is_on_y_identification_2 is_on_y_identification_2_object() const
  { return Is_on_y_identification_2(this, is_on_y_identification_op()); }

  Compare_x_at_limit_2 compare_x_at_limit_2_object() const
  { return Compare_x_at_limit_2(this, compare_x_at_limit_op()); }

  Compare_x_near_limit_2 compare_x_near_limit_2_object() const
  { return Compare_x_near_limit_2(this, compare_x_near_limit_op()); }

  Compare_x_on_boundary_2 compare_x_on_boundary_2_object() const
  { return Compare_x_on_boundary_2(this, compare_x_on_boundary_op()); }

  Compare_x_near_boundary_2 compare_x_near_boundary_2_object() const
  { return Compare_x_near_boundary_2(this, compare_x_near_boundary_op()); }

  //@}
};

template <class OutputStream>
OutputStream & operator<<(OutputStream & os, Comparison_result cr)
{
  os << ((cr == SMALLER) ? "SMALLER" : (cr == EQUAL) ? "EQUAL" : "LARGER");
  return os;
}

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
