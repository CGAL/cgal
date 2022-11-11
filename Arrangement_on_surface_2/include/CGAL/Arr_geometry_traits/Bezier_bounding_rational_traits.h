// Copyright (c) 2006,2007,2008,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Iddo Hanniel <iddoh@cs.technion.ac.il>
//                 Ron Wein     <wein@post.tau.ac.il>

#ifndef CGAL_BEZIER_BOUNDING_RATIONAL_TRAITS_H
#define CGAL_BEZIER_BOUNDING_RATIONAL_TRAITS_H

#include <CGAL/license/Arrangement_on_surface_2.h>


/*! \file
 * Definition of the Bezier_bounding_rational_traits<Kernel> class.
 */

#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Arr_geometry_traits/de_Casteljau_2.h>

#include <deque>
#include <list>
#include <math.h>

namespace CGAL {

/*! \struct _Bez_point_bound
 * Representation of a bounding interval for a point on a Bezier curve.
 * Basically, we store an interval [t_min, t_max] such that the point p
 * equals B(t) for some point t_min <= t <= t_max. We also store a bounding
 * polynomial for the point (obtain by de Casteljau's algorithm).
 */
template <class Kernel_>
struct _Bez_point_bound
{
  /*! \enum Type
   * The point type.
   */
  enum Type
  {
    RATIONAL_PT,                /*!< A point with rational coordinates .*/
    VERTICAL_TANGENCY_PT,       /*!< A vertical tangency point .*/
    INTERSECTION_PT,            /*!< An intersection point .*/
    UNDEFINED
  };

  typedef Kernel_                       Kernel;
  typedef typename Kernel::FT           NT;
  typedef typename Kernel::Point_2      Point_2;
  typedef std::deque<Point_2>           Control_points;

  Type              type;         /*!< The point type. */
  Control_points    ctrl;         /*!< The control point whose convex
                                       hull contains the point. */
  NT                t_min;        /*!< Minimal t value. */
  NT                t_max;        /*!< Maximal t value. */
  bool              can_refine;   /*!< Can we refine the current
                                       representation. */

  /*! Default constructor. */
  _Bez_point_bound() :
    type (UNDEFINED),
    ctrl(),
    t_min(), t_max(),
    can_refine (false)
  {}

  /*! Constructor with parameters. */
  _Bez_point_bound (Type _type,
                    const Control_points& _ctrl,
                    const NT& _t_min, const NT& _t_max,
                    bool _can_refine) :
    type (_type),
    ctrl (_ctrl),
    t_min (_t_min), t_max (_t_max),
    can_refine (_can_refine)
  {}
};

/*! \struct _Bez_point_bbox
 * A bounding box for a point on a Bezier curve.
 */
template <class Kernel_>
struct _Bez_point_bbox
{
  typedef Kernel_                       Kernel;
  typedef typename Kernel::FT           NT;
  typedef _Bez_point_bbox<Kernel>       Self;

  NT    min_x;              /*!< Lower bound for the x-coordinate. */
  NT    max_x;              /*!< Upper bound for the x-coordinate. */
  NT    min_y;              /*!< Lower bound for the y-coordinate. */
  NT    max_y;              /*!< Upper bound for the y-coordinate. */

  /*! Default constructor. */
  _Bez_point_bbox() :
    min_x(0), max_x(0),
    min_y(0), max_y(0)
  {}

  /*! Constructor with parameters. */
  _Bez_point_bbox (const NT& _min_x, const NT& _max_x,
                   const NT& _min_y, const NT& _max_y) :
    min_x(_min_x), max_x(_max_x),
    min_y(_min_y), max_y(_max_y)
  {}

  /*! Add two bounding boxes, and obtain a box bounding them both. */
  Self operator+ (const Self& other) const
  {
    return (Self ((std::min)(min_x, other.min_x), (std::max)(max_x, other.max_x),
                  (std::min)(min_y, other.min_y), (std::max)(max_y, other.max_y)));
  }

  /*! Addition and assignment. */
  void operator+= (const Self& other)
  {
    min_x = (std::min)(min_x, other.min_x);
    max_x = (std::max)(max_x, other.max_x);
    min_y = (std::min)(min_y, other.min_y);
    max_y = (std::max)(max_y, other.max_y);
    return;
  }

  /*! Check whether the two bounding boxed overlap. */
  bool overlaps (const Self& other) const
  {
    if ((CGAL::compare (max_x, other.min_x) == CGAL::SMALLER) ||
        (CGAL::compare (min_x, other.max_x) == CGAL::LARGER)  ||
        (CGAL::compare (max_y, other.min_y) == CGAL::SMALLER) ||
        (CGAL::compare (min_y, other.max_y) == CGAL::LARGER))
    {
      return (false);
    }
    return (true);
  }

  /*! Check whether the two bounding boxed overlap in their x-range. */
  bool overlaps_x (const Self& other) const
  {
    if ((CGAL::compare (max_x, other.min_x) == CGAL::SMALLER) ||
        (CGAL::compare (min_x, other.max_x) == CGAL::LARGER))
    {
      return (false);
    }
    return (true);
  }
};

/*!\ class Bezier_bounding_rational_traits
 * A traits class for performing bounding operations on Bezier curves,
 * assuming the NT is rational.
 */
template <typename _Kernel>
class Bezier_bounding_rational_traits
{
public:

  typedef _Kernel                                   Kernel;
  typedef typename Kernel::FT                       NT;

  typedef _Bez_point_bound<Kernel>                  Bez_point_bound;
  typedef _Bez_point_bbox<Kernel>                   Bez_point_bbox;
  typedef typename Bez_point_bound::Control_points  Control_points;

  /*! \struct Vertical_tangency_point
   * Representation of an approximated vertical tangency point.
   */
  struct Vertical_tangency_point
  {
    Bez_point_bound     bound;  /*!< A point-on-curve bound. */
    Bez_point_bbox      bbox;   /*!< Bounding box for the point. */

    /*! Default constructor. */
    Vertical_tangency_point () :
      bound (),
      bbox ()
    {}

    /*! Constructor. */
    Vertical_tangency_point (const Bez_point_bound& _bound,
                             const Bez_point_bbox& _bbox) :
      bound (_bound),
      bbox (_bbox)
    {}
  };

  /*! \struct Vertical_tangency_point
   * Representation of an approximated intersection point.
   */
  struct Intersection_point
  {
    Bez_point_bound     bound1; /*!< A point-on-curve bound for the
                                     first curve. */
    Bez_point_bound     bound2; /*!< A point-on-curve bound for the
                                     second curve. */
    Bez_point_bbox      bbox;   /*!< Bounding box for the point. */

    /*! Default constructor. */
    Intersection_point() :
      bound1(),
      bound2(),
      bbox()
    {}

    /*! Constructor. */
    Intersection_point (const Bez_point_bound& _bound1,
                        const Bez_point_bound& _bound2,
                        const Bez_point_bbox& _bbox) :
      bound1 (_bound1),
      bound2 (_bound2),
      bbox (_bbox)
    {}

  };

protected:

  typedef typename Kernel::Point_2                  Point_2;
  typedef typename Kernel::Vector_2                 Vector_2;
  typedef typename Kernel::Line_2                   Line_2;
  typedef typename Kernel::Direction_2              Direction_2;

  // Data members:
  NT            m_accuracy_bound;       /*!< An accuracy bound. */
  Kernel        m_kernel;               /*!< A geometry kernel. */
  unsigned int  m_active_nodes;         /*!< Limit for the number of recursive
                                             calls for the approximated
                                             intersection procedure. */

  // Kernel functors:
  typename Kernel::Construct_vector_2             f_construct_vector;
  typename Kernel::Construct_line_2               f_construct_line;
  typename Kernel::Construct_direction_2          f_construct_direction;
  typename Kernel::Construct_opposite_direction_2 f_opposite_direction;
  typename Kernel::Construct_translated_point_2   f_construct_point;
  typename Kernel::Equal_2                        f_equal;
  typename Kernel::Compare_x_2                    f_compare_x;
  typename Kernel::Compare_y_2                    f_compare_y;
  typename Kernel::Orientation_2                  f_orientation;
  typename Kernel::Counterclockwise_in_between_2  f_ccw_in_between;
  typename Kernel::Less_signed_distance_to_line_2 f_less_signed_distance;
  typename Kernel::Oriented_side_2                f_oriented_side;
  typename Kernel::Intersect_2                    f_intersect;

public:

  /*!
   * Constructor.
   * \param bound Accuracy bound. We cannot determine the order of two values
   *              whose absolute difference is smaller than this bound.
   */
  Bezier_bounding_rational_traits (double bound = 0.00000000000000011) :
    m_accuracy_bound (bound),
    m_active_nodes (0)
  {
    // Construct the kernel functors.
    f_construct_vector = m_kernel.construct_vector_2_object();
    f_construct_line = m_kernel.construct_line_2_object();
    f_construct_direction = m_kernel.construct_direction_2_object();
    f_opposite_direction = m_kernel.construct_opposite_direction_2_object();
    f_construct_point = m_kernel.construct_translated_point_2_object();
    f_equal = m_kernel.equal_2_object();
    f_compare_x = m_kernel.compare_x_2_object();
    f_compare_y = m_kernel.compare_y_2_object();
    f_orientation = m_kernel.orientation_2_object();
    f_ccw_in_between = m_kernel.counterclockwise_in_between_2_object();
    f_less_signed_distance = m_kernel.less_signed_distance_to_line_2_object();
    f_oriented_side = m_kernel.oriented_side_2_object();
    f_intersect = m_kernel.intersect_2_object();
  }

  /// \name Handling intersection points.
  //@{

  /*!
   * Check whether the curve specified by the given control polygon may
   * have self intersections.
   * \param cp The control polygon of the curve.
   * \return Whether the curve may be self-intersecting (but not necessarily).
   *         In case the function returns false, there can be absolutely no
   *         self-intersections.
   */
  bool may_have_self_intersections (const Control_points& cp)
  {
    // In case the control polygon is convex, the Bezier curve cannot be
    // self-intersecting.
    if (is_convex_2 (cp.begin(), cp.end(), m_kernel))
      return (false);

    // In case the angular span of the control polygon is less than 180
    // degrees, it cannot be self-intersecting.
    Vector_2    v_min;
    Vector_2    v_max;

    if (_compute_angular_span (cp, v_min, v_max))
      return (false);

    // Otherwise, the curve may be self-intersecting.
    return (true);
  }

  /*!
   * Compute bounds for the intersection points between two Bezier curves.
   * \param cp1 The control polygon of the first curve.
   * \param cp2 The control polygon of the second curve.
   * \param oi Output: The bounds of the intersection points.
   *                   The value-type of this iterator is Intersection_point.
   */
  template <class OutputIterator>
  OutputIterator compute_intersection_points (const Control_points& cp1,
                                              const Control_points& cp2,
                                              OutputIterator oi)
  {
    // Call the recursive function.
    std::set<NT>                    dummy;
    std::list<Intersection_point>   ipts;

    m_active_nodes = 1;
    _compute_intersection_points (cp1, 0, 1,
                                  cp2, 0, 1,
                                  true, dummy,      // Check span.
                                  ipts);

    // Copy the computed points to the output iterator.
    typename std::list<Intersection_point>::iterator     it;

    for (it = ipts.begin(); it != ipts.end(); ++it)
    {
        *oi = *it;
        ++oi;
    }

    return (oi);
  }

  /*!
   * Refine the approximation of a given intersection point.
   * \param in_pt The intersection point.
   * \param ref_pt Output: The refined representation of the point.
   * \return Whether the point has been successfully refined (if not, the
   *         output point is the same as the input).
   */
  bool refine_intersection_point (const Intersection_point& in_pt,
                                  Intersection_point& ref_pt)
  {
    // Check if it is not possible to refine.
    if (! in_pt.bound1.can_refine || ! in_pt.bound2.can_refine)
    {
      ref_pt = in_pt;
      return (false);
    }

    // In case the point is already rational, there is not point in refining
    // it.
    if (in_pt.bound1.type == Bez_point_bound::RATIONAL_PT ||
        in_pt.bound2.type == Bez_point_bound::RATIONAL_PT)
    {
      ref_pt = in_pt;
      return (false);
    }

    // Check whether we can refine the parametric ranges of the two
    // originating Bezier curves.
    const Control_points&  cp1 = in_pt.bound1.ctrl;
    const NT&              t_min1 = in_pt.bound1.t_min;
    const NT&              t_max1 = in_pt.bound1.t_max;
    bool                   can_refine1 = can_refine (cp1, t_min1, t_max1);

    const Control_points&  cp2 = in_pt.bound2.ctrl;
    const NT&              t_min2 = in_pt.bound2.t_min;
    const NT&              t_max2 = in_pt.bound2.t_max;
    bool                   can_refine2 = can_refine (cp2, t_min2, t_max2);

    if (! can_refine1 || ! can_refine2)
    {
      // Failed in the bounded approach - stop the subdivision and indicate
      // that the subdivision has failed.
      ref_pt = in_pt;
      ref_pt.bound1.can_refine = false;
      ref_pt.bound2.can_refine = false;

      return (false);
    }

    // Apply de Casteljau's algorithm and bisect both bounding polylines.
    Control_points         cp1a, cp1b;
    const NT               t_mid1 = (t_min1 + t_max1) / 2;

    bisect_control_polygon_2 (cp1.begin(), cp1.end(),
                              std::back_inserter(cp1a),
                              std::front_inserter(cp1b));

    Control_points         cp2a, cp2b;
    const NT               t_mid2 = (t_min2 + t_max2) / 2;

    bisect_control_polygon_2 (cp2.begin(), cp2.end(),
                              std::back_inserter(cp2a),
                              std::front_inserter(cp2b));

    // Catch situations that we cannot further refine, in order to prevent
    // having multiple points as the refined point.
    const bool can_refine_pt_cp1a = can_refine (cp1a, t_min1, t_mid1);
    const bool can_refine_pt_cp1b = can_refine (cp1b, t_mid1, t_max1);
    const bool can_refine_pt_cp2a = can_refine (cp2a, t_min2, t_mid2);
    const bool can_refine_pt_cp2b = can_refine (cp2b, t_mid2, t_max2);

    if (! can_refine_pt_cp1a || ! can_refine_pt_cp1b ||
        ! can_refine_pt_cp2a || ! can_refine_pt_cp2b)
    {
      // Construct an inconclusive point bound, which includes all parameters
      // found so far, with can_refine = false, and with the bounding
      // box of cp1.
      Bez_point_bound   bound1 (Bez_point_bound::INTERSECTION_PT,
                                cp1, t_min1, t_max1,
                                false);         // Cannot refine further.
      Bez_point_bound   bound2 (Bez_point_bound::INTERSECTION_PT,
                                cp2, t_min2, t_max2,
                                false);         // Cannot refine further.
      Bez_point_bbox    ref_bbox;

      construct_bbox (cp1, ref_bbox);
      ref_pt = Intersection_point (bound1, bound2, ref_bbox);
      return (true);
    }

    // Use the bisection method to refine the intersection point.
    // We assume that there is no need to check the span, as the input point
    // already has all the necessary data.
    std::list<Intersection_point>  ipts;
    std::set<NT>                   dummy;

    m_active_nodes = 1;
    _compute_intersection_points (cp1a, t_min1, t_mid1,
                                  cp2a, t_min2, t_mid2,
                                  false, dummy,     // Don't check span.
                                  ipts);
    if (! ipts.empty())
    {
      if (ipts.size() == 1)
      {
        ref_pt = ipts.front();
        return (true);
      }
      ipts.clear();
    }

    m_active_nodes = 1;
    _compute_intersection_points (cp1a, t_min1, t_mid1,
                                  cp2b, t_mid2, t_max2,
                                  false, dummy,     // Don't check span.
                                  ipts);
    if (! ipts.empty())
    {
      if (ipts.size() == 1)
      {
        ref_pt = ipts.front();
        return (true);
      }
      ipts.clear();
    }

    m_active_nodes = 1;
    _compute_intersection_points (cp1b, t_mid1, t_max1,
                                  cp2a, t_min2, t_mid2,
                                  false, dummy,     // Don't check span.
                                  ipts);
    if (! ipts.empty())
    {
      if (ipts.size() == 1)
      {
        ref_pt = ipts.front();
        return (true);
      }
      ipts.clear();
    }

    m_active_nodes = 1;
    _compute_intersection_points (cp1b, t_mid1, t_max1,
                                  cp2b, t_mid2, t_max2,
                                  false, dummy,     // Don't check span.
                                  ipts);

    CGAL_assertion (ipts.size() == 1);
    ref_pt = ipts.front();
    return (true);
  }
  //@}

  /// \name Handling vertcial tangency points.
  //@{

  /*!
   * Compute bounds for the vertical tangency points of a Bezier curves.
   * The points are returned sorted by their parameter value.
   * \param cp The control polygon of the curve.
   * \param oi Output: The bounds of the vertical tangency points.
   *                   The value-type is Vertical_tangency_point.
   */
  template <class OutputIterator>
  OutputIterator
      compute_vertical_tangency_points (const Control_points& cp,
                                        OutputIterator oi)
  {
    // Call the recursive function on the entire curve (0 <= t <= 1).
    std::list<Vertical_tangency_point>            vpts;

    _compute_vertical_tangency_points (cp,
                                       0, 1,
                                       vpts);

    // Copy the computed points to the output iterator, sorted by their
    // parameter value (note we compare the t_min values of the bounds; this
    // is possible as different vertical tangency points are expected to
    // have disjoint parametric ranges).
    typename std::list<Vertical_tangency_point>::iterator  vpt_it;
    typename std::list<Vertical_tangency_point>::iterator  vpt_end;
    typename std::list<Vertical_tangency_point>::iterator  vpt_min;

    while (! vpts.empty())
    {
      // Locate the vertical tangency point with minimal t-value.
      vpt_min = vpt_it = vpts.begin();
      vpt_end = vpts.end();

      ++vpt_it;
      while (vpt_it != vpt_end)
      {
        if (CGAL::compare (vpt_it->bound.t_min,
                           vpt_min->bound.t_min) == SMALLER)
        {
          vpt_min = vpt_it;
        }

        ++vpt_it;
      }

      // Copy this point to the output iterator.
      *oi = *vpt_min;
      ++oi;

      // Remove it from the list.
      vpts.erase (vpt_min);
    }

    return (oi);
  }

  /*!
   * Refine the approximation of a given vertical tangency point.
   * \param in_pt The vertical tangency point.
   * \param ref_pt Output: The refined representation of the point.
   * \return Whether the point has been successfully refined (if not, the
   *         output point is the same as the input).
   */
  bool refine_vertical_tangency_point (const Vertical_tangency_point& in_pt,
                                       Vertical_tangency_point& ref_pt)
  {
    // Check if it is not possible to refine.
    if (! in_pt.bound.can_refine ||
        in_pt.bound.type == Bez_point_bound::RATIONAL_PT)
    {
      ref_pt = in_pt;
      return (false);
    }

    // Check whether we can refine the parametric ranges of the originating
    // Bezier curve.
    const Control_points&  cp = in_pt.bound.ctrl;
    const NT&              t_min = in_pt.bound.t_min;
    const NT&              t_max = in_pt.bound.t_max;

    if (! can_refine (cp, t_min, t_max))
    {
      // Failed in the bounded approach - stop the subdivision and indicate
      // that the subdivision has failed.
      ref_pt = in_pt;
      ref_pt.bound.can_refine = false;

      return (false);
    }

    // Apply de Casteljau's algorithm and bisect the bounding polyline.
    Control_points         cp_a, cp_b;
    const NT               t_mid = (t_min + t_max) / 2;

    bisect_control_polygon_2 (cp.begin(), cp.end(),
                              std::back_inserter(cp_a),
                              std::front_inserter(cp_b));

    // Handle the case where t_mid is a vertical (rational) tangency
    // point, by checking whether the first vector of right is vertical.
    if (f_compare_x (cp_b[0], cp_b[1]) == EQUAL)
    {
      // The subdivision point at t_mid is a vertical tangency point with
      // rational coordinates.
      const typename Kernel::Point_2  &vpt = cp_b[0];

      Bez_point_bound   bound (Bez_point_bound::RATIONAL_PT,
                               cp_b, t_mid, t_mid, true);
      Bez_point_bbox    bbox (vpt.x(), vpt.x(), vpt.y(), vpt.y());

      ref_pt = Vertical_tangency_point (bound, bbox);
      return (true);
    }

    const bool can_refine_pt_cp_a = can_refine (cp_a, t_min, t_mid);
    const bool can_refine_pt_cp_b = can_refine (cp_b, t_mid, t_max);

    if (! can_refine_pt_cp_a || ! can_refine_pt_cp_b)
    {
      // Construct an inconclusive point bound, which includes all parameters
      // found so far, with can_refine = false, and with the bounding
      // box of the curve.
      Bez_point_bound   bound (Bez_point_bound::VERTICAL_TANGENCY_PT,
                               cp, t_min, t_max,
                               false);         // Cannot refine further.
      Bez_point_bbox    ref_bbox;

      construct_bbox (cp, ref_bbox);
      ref_pt = Vertical_tangency_point (bound, ref_bbox);
      return (true);
    }

    // Compute the vertical tangency points of the two subcurves in order to
    // refine the vertical tangency point.
    std::list<Vertical_tangency_point>  vpts;

    _compute_vertical_tangency_points(cp_a, t_min, t_mid, vpts);
    _compute_vertical_tangency_points(cp_b, t_mid, t_max, vpts);

    CGAL_assertion(vpts.size() == 1);

    ref_pt = vpts.front();
    return (true);
  }
  //@}

  /// \name Predicates.
  //@{

  /*!
   * Check whether the current polyline can be refined.
   * In this case we simply check the values that define the parametric
   * range are well-separated.
   */
  bool can_refine (const Control_points& ,
                   const NT& t_min, const NT& t_max)
  {
    return (CGAL::compare (t_max - t_min, m_accuracy_bound) != SMALLER);
  }

  /*!
   * Compare the slopes of two Bezier curve at their intersection point,
   * given the point-on-curve bounds at this point.
   * \param bound1 The point-on-curve bound of the first curve.
   * \param bound2 The point-on-curve bound of the second curve.
   * \return The comparison result.
   */
  Comparison_result
  compare_slopes_at_intersection_point (const Bez_point_bound& bound1,
                                        const Bez_point_bound& bound2)
  {
    const Control_points& cp1 = bound1.ctrl;
    const Control_points& cp2 = bound2.ctrl;

    // An (expensive) check that the angular spans do not overlap.
    CGAL_expensive_precondition_code (
      Vector_2    v_min1;
      Vector_2    v_max1;
      const bool  span_ok1 = _compute_angular_span (cp1,
                                                    v_min1, v_max1);
      Vector_2    v_min2;
      Vector_2    v_max2;
      const bool  span_ok2 = _compute_angular_span (cp2,
                                                    v_min2, v_max2);
      bool        spans_overlap = true;

      if (span_ok1 && span_ok2)
      {
          spans_overlap = _angular_spans_overlap (v_min1, v_max1,
                                                  v_min2, v_max2);
      }
    );

    CGAL_expensive_precondition (! spans_overlap);

    // As the angular spans of the control polygons do not overlap, we can
    // just compare any vector of the hodograph spans.
    const Vector_2  dir1 = f_construct_vector (cp1.front(), cp1.back());
    const Vector_2  dir2 = f_construct_vector (cp2.front(), cp2.back());

    // The slopes are given by dir1.y() / dir1.x() and dir2.y()/dir2.x().
    // However, to avoid potential division by zero, we compare:
    // dir1.y()*dir2.x() and dir2.y()*dir1.x(), after considering the signs
    // of the xs.
    Comparison_result  res;
    bool               swap_res = false;

    if (CGAL::sign (dir1.x()) == CGAL::NEGATIVE)
        swap_res = ! swap_res;

    if (CGAL::sign (dir2.x()) == CGAL::NEGATIVE)
        swap_res = ! swap_res;

    res = CGAL::compare (dir1.y() * dir2.x(), dir2.y() * dir1.x());
    if (swap_res)
        res = CGAL::opposite (res);

    return (res);
  }

  /*!
   * Construct a bounding box for the given control polygon.
   * \param cp A sequence of control point (the control polgon).
   * \param bbox Output: The bounding box.
   * \pre cp is not empty.
   */
  void construct_bbox (const Control_points& cp,
                       Bez_point_bbox& bez_bbox)
  {
    CGAL_precondition (! cp.empty());

    // Go over the points, and locate the ones with extremal x and y values.
    typename Control_points::const_iterator   it = cp.begin();
    typename Control_points::const_iterator   it_end = cp.end();
    typename Control_points::const_iterator   min_x = it;
    typename Control_points::const_iterator   max_x = it;
    typename Control_points::const_iterator   min_y = it;
    typename Control_points::const_iterator   max_y = it;

    while (it != it_end)
    {
      if (f_compare_x (*it, *min_x) == SMALLER)
        min_x = it;
      else if (f_compare_x (*it, *max_x) == LARGER)
        max_x = it;

      if (f_compare_y (*it, *min_y) == SMALLER)
        min_y = it;
      else if (f_compare_y (*it, *max_y) == LARGER)
        max_y = it;

      ++it;
    }

    // Set the bounding box.
    bez_bbox.min_x = min_x->x();
    bez_bbox.max_x = max_x->x();
    bez_bbox.min_y = min_y->y();
    bez_bbox.max_y = max_y->y();
    return;
  }
  //@}

private:

  /*!
   * Check whether the given control polygon defines an x-monotone (or
   * y-monotone) Bezier curve.
   * \param cp The control points.
   * \param check_x (true) to check x-monotonicity;
   *                (false) to check y-monotonicity.
   * \return Is the curve x-monotone (or y-monotone).
   */
  bool _is_monotone (const Control_points& cp, bool check_x) const
  {
    Comparison_result                  curr_res;
    Comparison_result                  res = EQUAL;

    // Look for the first pair of consecutive points whose x-coordinate
    // (or y-coordinate) are not equal. Their comparsion result will be
    // set as the "reference" comparison result.
    typename Control_points::const_iterator  pt_curr = cp.begin();
    typename Control_points::const_iterator  pt_end = cp.end();
    typename Control_points::const_iterator  pt_next = pt_curr;

    ++pt_next;
    while (pt_next != pt_end)
    {
      curr_res = check_x ? f_compare_x (*pt_curr, *pt_next) :
                           f_compare_y (*pt_curr, *pt_next);

      pt_curr = pt_next;
      ++pt_next;

      if (curr_res != EQUAL)
      {
        res = curr_res;
        break;
      }
    }

    if (pt_next == pt_end || res == EQUAL)
    {
      // If we reached here, we have an iso-parallel segment.
      // In our context, we consider a horizontal segment as non-y-monotone
      // but a vertical segment is considered x-monotone.
      if (check_x)
        return (true);  // Vertical segment when checking x-monoticity

      return (false);  // Horizontal segment when checking y-monoticity
    }

    // Go over the rest of the control polygons, and make sure that the
    // comparison result between each pair of consecutive points is the
    // same as the reference result.
    while (pt_next != pt_end)
    {
      curr_res = check_x ? f_compare_x (*pt_curr, *pt_next) :
                           f_compare_y (*pt_curr, *pt_next);

      if (curr_res != EQUAL && res != curr_res)
        return (false);

      pt_curr = pt_next;
      ++pt_next;
    }

    // If we reached here, the control polygon is x-monotone (or y-monotone).
    return (true);
  }

  /*!
   * Check whether the  control polygon defines an x-monotone Bezier curve.
   * \param cp The control points.
   * \return Is the curve x-monotone.
   */
  inline bool _is_x_monotone (const Control_points& cp) const
  {
    return (_is_monotone (cp, true));
  }

  /*!
   * Check whether the control polygon defines a y-monotone Bezier curve.
   * \param cp The control points.
   * \return Is the curve y-monotone.
   */
  inline bool _is_y_monotone (const Control_points& cp) const
  {
    return (_is_monotone (cp, false));
  }

  /*!
   * Compute the angular span of a Bezier curve, given by its control polygon.
   * \param cp The control points.
   * \param v_min Output: The minimal direction of the angular span.
   * \param v_max Output: The maximal direction of the angular span.
   * \return Whether the span is less than 90 degrees
   *         (as we don't want to work with larger spans).
   */
  bool _compute_angular_span (const Control_points& cp,
                              Vector_2& v_min,
                              Vector_2& v_max)
  {
    // Initialize the output directions.
    Vector_2    dir = f_construct_vector (cp.front(), cp.back());
    Vector_2    v;

    v_min = v_max = dir;

    // Go over the control points and examine the vectors defined by pairs
    // of consecutive control points.
    typename Control_points::const_iterator  pt_curr = cp.begin();
    typename Control_points::const_iterator  pt_end = cp.end();
    typename Control_points::const_iterator  pt_next = pt_curr;

    ++pt_next;
    while (pt_next != pt_end)
    {
      v = f_construct_vector (*pt_curr, *pt_next);

      // If the current vector forms a right-turn with dir, it is a candidate
      // for the minimal vector, otherwise it is a candidate for the maximal
      // vector in the span.
      bool                                   updated_span = false;

      if (f_orientation (f_construct_point (ORIGIN, v),
                         ORIGIN,
                         f_construct_point (ORIGIN, dir)) == RIGHT_TURN)
      {
        // Check if we need to update the minimal vector in the span.
        if (f_orientation (f_construct_point (ORIGIN, v),
                           ORIGIN,
                           f_construct_point (ORIGIN, v_min)) == RIGHT_TURN)
        {
          v_min = v;
          updated_span = true;
        }
      }
      else
      {
        // Check if we need to update the maximal vector in the span.
        if (f_orientation (f_construct_point (ORIGIN, v),
                           ORIGIN,
                           f_construct_point (ORIGIN, v_max)) == LEFT_TURN)
        {
          v_max = v;
          updated_span = true;
        }
      }

      // Before proceeding, check if the current span is larger than
      // 180 degrees. We do that by checking that (v_min, 0, v_max)
      // is still a right-turn, and did not become a left-turn.
      if (updated_span &&
          f_orientation (f_construct_point (ORIGIN, v_min),
                         ORIGIN,
                         f_construct_point (ORIGIN, v_max)) != RIGHT_TURN)
      {
        return (false);
      }

      // Move to the next pair of points.
      pt_curr = pt_next;
      ++pt_next;
    }

    // If we reached here, the angular span is less than 180 degrees.
    return (true);
  }

  /*!
   * Check whether to angular spans overlap.
   * \param v_min1 The minimal vector in the first angular span.
   * \param v_max1 The maximal vector in the first angular span.
   * \param v_min2 The minimal vector in the second angular span.
   * \param v_max2 The maximal vector in the second angular span.
   * \return Do the two spans overlap.
   */
  bool _angular_spans_overlap (const Vector_2& v_min1,
                               const Vector_2& v_max1,
                               const Vector_2& v_min2,
                               const Vector_2& v_max2)
  {
    const Direction_2     dir_l1 = f_construct_direction (v_min1);
    const Direction_2     dir_r1 = f_construct_direction (v_max1);
    const Direction_2     dir_l2 = f_construct_direction (v_min2);
    const Direction_2     dir_r2 = f_construct_direction (v_max2);

    //use the directions to test equality of spans
    if ( f_equal (dir_l1, dir_l2) && f_equal (dir_r1, dir_r2)) //spans are identical
      return (true);

    if ( f_equal (dir_l1, dir_r1) && f_equal (dir_r2, dir_l2) && f_equal (dir_l1, -dir_l2)) //spans are identical but opposite
      return (true);

    // Check whether any of the vectors of the second span (or their
    // opposite vectors) is between the vectors of the first span.
    if (! f_equal (dir_l1, dir_r1) &&
        (f_ccw_in_between (dir_l2, dir_l1, dir_r1) ||
         f_ccw_in_between (f_opposite_direction (dir_l2), dir_l1, dir_r1) ||
         f_ccw_in_between (dir_r2, dir_l1, dir_r1) ||
         f_ccw_in_between (f_opposite_direction (dir_r2), dir_l1, dir_r1)))
    {
      return (true);
    }

    // Check whether the left vector of the first span (or its opposite
    // vector) is between the vectors of the second span. Note that at this
    // point, the only possibility is that the first span is totally contained
    // within the second, so we do not need to check both vectors.
    if (! f_equal (dir_l2, dir_r2) &&
        (f_ccw_in_between (dir_l1, dir_l2, dir_r2) ||
         f_ccw_in_between (f_opposite_direction (dir_l1), dir_l2, dir_r2)))
    {
      return (true);
    }

    // If we reached here, the two angular spans do not overlap.
    return (false);
  }

  /*!
   * Construct a skewed bounding box for the control polygon.
   * The skewed bounding box is represented as a pair of lines, both parallel
   * to the line connecting the first and last points of the control polygon,
   * that bound the curve from above and below.
   * \param cp The control points.
   * \param l_min Output: The line with minimal signed distance.
   * \param l_max Output: The line with maximal signed distance.
   */
  void _skewed_bbox (const Control_points& cp,
                     Line_2& l_min, Line_2& l_max)
  {
    // Construct the line that connects the first and the last control points.
    const Line_2    l = f_construct_line (cp.front(), cp.back());
    const Vector_2  v = f_construct_vector (cp.front(), cp.back());

    // Select the points with minimum and maximum distance from the line l,
    // and construct two lines parallel to l that pass through these points.
    typename Control_points::const_iterator   pt_curr = cp.begin();
    typename Control_points::const_iterator   pt_end = cp.end();

    typename Control_points::const_iterator   pt_min = pt_curr;
    typename Control_points::const_iterator   pt_max = pt_curr;

    ++pt_curr;
    while (pt_curr != pt_end)
    {
      if (f_less_signed_distance (l, *pt_curr, *pt_min))
      {
        pt_min = pt_curr;
      }
      else if (f_less_signed_distance (l, *pt_max, *pt_curr))
      {
        pt_max = pt_curr;
      }

      ++pt_curr;
    }

    // Construct the output lines.
    l_min = f_construct_line (*pt_min, v);
    l_max = f_construct_line (*pt_max, v);
    return;
  }

  /*!
   * Check whether the endpoints of the two curves coincide.
   * If so, it adds the Intersection_point to the list if the endpoint hasn't
   * already been inserted there.
   * Note that this function is not able to detect an endpoint that lies in
   * the interior of the other curve. Such cases will not be approximated
   * and we will have to resort to exact computation.
   * \param cp1 The control points of the first curve.
   * \param t_min1 The lower bound of the parameter range of the first curve.
   * \param t_max1 The upper bound of the parameter range of the first curve.
   * \param cp2 The control points of the second curve.
   * \param t_min2 The lower bound of the parameter range of the second curve.
   * \param t_max2 The upper bound of the parameter range of the second curve.
   * \param iept_params Input/Output: The set of parameter values (for cp1)
   *                                  of known intersections at the endpoints.
   * \param ipts Input/Output: The computed intersection points.
   * \return Whether a pair of curve endpoints coincide.
   */
  bool _endpoints_coincide (const Control_points& cp1,
                            const NT& t_min1, const NT& t_max1,
                            const Control_points& cp2,
                            const NT& t_min2, const NT& t_max2,
                            std::set<NT>& iept_params,
                            std::list<Intersection_point>& ipts) const
  {
    // Get the first and last control points of each curve.
    const Point_2&  s1 = cp1.front();
    const Point_2&  t1 = cp1.back();
    const Point_2&  s2 = cp2.front();
    const Point_2&  t2 = cp2.back();

    // Check whether any pair of these endpoints conincide.
    NT      x, y;               // Coordinate of a common endpoint.
    NT      t_val1, t_val2;     // Its respective parameters.

    if (f_equal (s1, s2))
    {
      x = s1.x();
      y = s1.y();
      t_val1 = t_min1;
      t_val2 = t_min2;
    }
    else if (f_equal (s1, t2))
    {
      x = s1.x();
      y = s1.y();
      t_val1 = t_min1;
      t_val2 = t_max2;
    }
    else if (f_equal (t1, s2))
    {
      x = t1.x();
      y = t1.y();
      t_val1 = t_max1;
      t_val2 = t_min2;
    }
    else if (f_equal (t1, t2))
    {
      x = t1.x();
      y = t1.y();
      t_val1 = t_max1;
      t_val2 = t_max2;
    }
    else
    {
      // No common endpoint found:
      return (false);
    }

    // Try inserting the parameter value t1 into the set of parameters of
    // already discovered intersecting endpoints.
    std::pair<typename std::set<NT>::iterator,
              bool>                         res = iept_params.insert (t_val1);

    if (res.second)
    {
      // In case the insertion has succeeded, report on a new rational
      // intersection point.
      Bez_point_bound   bound1 (Bez_point_bound::RATIONAL_PT,
                                cp1, t_val1, t_val1, true);
      Bez_point_bound   bound2 (Bez_point_bound::RATIONAL_PT,
                                cp2, t_val2, t_val2, true);
      Bez_point_bbox    bbox (x, x, y, y);

      ipts.push_back (Intersection_point (bound1, bound2, bbox));
    }

    return (true);
  }

  /*!
   * An auxilary recursive function for computing the approximated
   * intersection points between two Bezier curves.
   * \param cp1 The control points of the first curve.
   * \param t_min1 The lower bound of the parameter range of the first curve.
   * \param t_max1 The upper bound of the parameter range of the first curve.
   * \param cp2 The control points of the second curve.
   * \param t_min2 The lower bound of the parameter range of the second curve.
   * \param t_max2 The upper bound of the parameter range of the second curve.
   * \param check_span Should we check the angular span of the curves.
   * \param iept_params Input/Output: The set of parameter values (for cp1)
   *                                  of known intersections at the endpoints.
   * \param ipts Input/Output: The computed intersection points.
   */
  void _compute_intersection_points (const Control_points& cp1,
                                     const NT& t_min1, const NT& t_max1,
                                     const Control_points& cp2,
                                     const NT& t_min2, const NT& t_max2,
                                     bool check_span,
                                     std::set<NT>& iept_params,
                                     std::list<Intersection_point>& ipts)
  {
    // Check if we got to subdivision termination criteria.
    const bool  can_refine1 = can_refine (cp1, t_min1, t_max1);
    const bool  can_refine2 = can_refine (cp2, t_min2, t_max2);;

    if (! can_refine1 || ! can_refine2)
    {
      // It is not possible to further refine the approximation, so we stop
      // the recursive subdivision, and construct an output intersection point
      // with can_refine = false. Note that we do not need a tight bounding
      // box here, since we cannot refine it anyway, so we take cp1's bbox
      // (which certainly contains the intersection point).
      Bez_point_bound   bound1 (Bez_point_bound::INTERSECTION_PT,
                                cp1, t_min1, t_max1,
                                false);         // Cannot refine further.
      Bez_point_bound   bound2 (Bez_point_bound::INTERSECTION_PT,
                                cp2, t_min2, t_max2,
                                false);         // Cannot refine further.
      Bez_point_bbox    ipt_bbox;

      construct_bbox (cp1, ipt_bbox);
      ipts.push_back (Intersection_point (bound1, bound2, ipt_bbox));

      m_active_nodes--;
      return;
    }

    // Construct bounding boxes for the two curves and check whether they
    // overlap.
    Bez_point_bbox      bbox1;
    Bez_point_bbox      bbox2;

    construct_bbox (cp1, bbox1);
    construct_bbox (cp2, bbox2);

    if (! bbox1.overlaps (bbox2))
    {
      // The bounding boxes do not overlap, so the two input curves do not
      // intersect.
      m_active_nodes--;
      return;
    }

    // Check the angular spans, if necessary.
    bool        spans_overlap = false;

    if (check_span)
    {
      // Compute the angular spans of the two curves.
      Vector_2    v_min1, v_max1;
      const bool  span_ok1 = _compute_angular_span (cp1,
                                                    v_min1, v_max1);
      Vector_2    v_min2, v_max2;
      const bool  span_ok2 = _compute_angular_span (cp2,
                                                    v_min2, v_max2);

      if (span_ok1 && span_ok2)
      {
        spans_overlap = _angular_spans_overlap (v_min1, v_max1,
                                                v_min2, v_max2);
      }
      else
      {
        // One of the spans is greater than 180 degrees, so we do not have to
        // check for overlaps (we know that it must overlap the other span).
        spans_overlap = true;
      }
    }

    if (! spans_overlap)
    {
      // In case the spans do not overlap, we potentially have a single
      // intersection point.

      // Checking for endpoint intersections.
      if (_endpoints_coincide (cp1, t_min1, t_max1,
                               cp2, t_min2, t_max2,
                               iept_params,
                               ipts))
      {
        // As we have located the single intersection point (a common endpoint
        // in this case), we can stop here.
        m_active_nodes--;
        return;
      }

      // Construct the skewed bounding boxes for the two curves and check
      // whether the endpoints of the first curve lie on opposite side of
      // the skewed bounding box of the second curve, and vice versa.
      Line_2               skew1a, skew1b;
      Line_2               skew2a, skew2b;
      const Point_2&       s1 = cp1.front();
      const Point_2&       t1 = cp1.back();
      const Point_2&       s2 = cp2.front();
      const Point_2&       t2 = cp2.back();

      _skewed_bbox(cp1, skew1a, skew1b);
      _skewed_bbox(cp2, skew2a, skew2b);

      const Oriented_side  or_2a_s1 = f_oriented_side (skew2a, s1);
      const Oriented_side  or_2a_t1 = f_oriented_side (skew2a, t1);

      const Oriented_side  or_2b_s1 = f_oriented_side (skew2b, s1);
      const Oriented_side  or_2b_t1 = f_oriented_side (skew2b, t1);

      const bool           s1_t1_are_opposite =
          ((or_2a_s1 == CGAL::opposite (or_2a_t1)) &&
           (or_2b_s1 == CGAL:: opposite (or_2b_t1))
          ) ||
          ( skew2a==skew2b && or_2a_s1!=or_2a_t1 );

      const Oriented_side  or_1a_s2 = f_oriented_side (skew1a, s2);
      const Oriented_side  or_1a_t2 = f_oriented_side (skew1a, t2);

      const Oriented_side  or_1b_s2 = f_oriented_side (skew1b, s2);
      const Oriented_side  or_1b_t2 = f_oriented_side (skew1b, t2);

      const bool           s2_t2_are_opposite =
          ((or_1a_s2 == CGAL::opposite (or_1a_t2)) &&
           (or_1b_s2 == CGAL::opposite (or_1b_t2)))||
          ( skew1a==skew1b && or_1a_s2!=or_1a_t2 );

      if (s1_t1_are_opposite && s2_t2_are_opposite)
      {
        // Construct a finer bounding box for the intersection point from
        // the intersection of the two skewed bounding boxes.
        Bez_point_bbox  ipt_bbox;
        Control_points  aux_vec;

        auto res1 = f_intersect(skew1a, skew2a);
        const Point_2* p1 = boost::get<Point_2>(&*res1);
        if (! p1) CGAL_error();
        aux_vec.push_back(*p1);

        auto res2 = f_intersect(skew1a, skew2b);
        const Point_2* p2 = boost::get<Point_2>(&*res2);
        if (! p2) CGAL_error();
        aux_vec.push_back(*p2);

        auto res3 = f_intersect(skew1b, skew2a);
        const Point_2* p3 = boost::get<Point_2>(&*res3);
        if (! p3) CGAL_error();
        aux_vec.push_back(*p3);

        auto res4 = f_intersect (skew1b, skew2b);
        const Point_2* p4 = boost::get<Point_2>(&*res4);
        if (! p4) CGAL_error();
        aux_vec.push_back(*p4);

        construct_bbox (aux_vec, ipt_bbox);

        // Report on the intersection point we have managed to isolate.
        Bez_point_bound     bound1 (Bez_point_bound::INTERSECTION_PT,
                                    cp1, t_min1, t_max1,
                                    true);      // We can further refine it.
        Bez_point_bound     bound2 (Bez_point_bound::INTERSECTION_PT,
                                    cp2, t_min2, t_max2,
                                    true);      // We can further refine it.

        ipts.push_back (Intersection_point (bound1, bound2, ipt_bbox));

        m_active_nodes--;
        return;
      }
    }

    // Apply de Casteljau's algorithm and bisect both bounding polylines.
    Control_points         cp1a, cp1b;
    const NT               t_mid1 = (t_min1 + t_max1) / 2;

    bisect_control_polygon_2 (cp1.begin(), cp1.end(),
                              std::back_inserter(cp1a),
                              std::front_inserter(cp1b));

    Control_points         cp2a, cp2b;
    const NT               t_mid2 = (t_min2 + t_max2) / 2;

    bisect_control_polygon_2 (cp2.begin(), cp2.end(),
                              std::back_inserter(cp2a),
                              std::front_inserter(cp2b));

    // Recursively compute the intersection points on all pairs of bisected
    // curves.
    m_active_nodes += 4;
    if (m_active_nodes > 4 * cp1.size() * cp2.size())
    {
      // It is not possible to further refine the approximation, as the
      // number of active nodes in the recursive search is too high. We stop
      // the recursive subdivision, and construct an output intersection point
      // with can_refine = false. Note that we do not need a tight bounding
      // box here, since we cannot refine it anyway, so we take cp1's bbox
      // (which certainly contains the intersection point).
      Bez_point_bound   bound1 (Bez_point_bound::INTERSECTION_PT,
                                cp1, t_min1, t_max1,
                                false);         // Cannot refine further.
      Bez_point_bound   bound2 (Bez_point_bound::INTERSECTION_PT,
                                cp2, t_min2, t_max2,
                                false);         // Cannot refine further.
      Bez_point_bbox    ipt_bbox;

      construct_bbox (cp1, ipt_bbox);
      ipts.push_back (Intersection_point (bound1, bound2, ipt_bbox));

      return;
    }

    _compute_intersection_points (cp1a, t_min1, t_mid1,
                                  cp2a, t_min2, t_mid2,
                                  spans_overlap,
                                  iept_params,
                                  ipts);

    _compute_intersection_points (cp1a, t_min1, t_mid1,
                                  cp2b, t_mid2, t_max2,
                                  spans_overlap,
                                  iept_params,
                                  ipts);

    _compute_intersection_points (cp1b, t_mid1, t_max1,
                                  cp2a, t_min2, t_mid2,
                                  spans_overlap,
                                  iept_params,
                                  ipts);

    _compute_intersection_points (cp1b, t_mid1, t_max1,
                                  cp2b, t_mid2, t_max2,
                                  spans_overlap,
                                  iept_params,
                                  ipts);

    m_active_nodes--;
    return;
  }

  /*!
   * An auxilary recursive function for computing the approximated vertical
   * tangency points of a Bezier curves.
   * \param cp The control points of the curve.
   * \param t_min The lower bound of the parameter range of the curve.
   * \param t_max The upper bound of the parameter range of the curve.
   * \param vpts Input/Output: The computed vertical tangency points.
   */
  void _compute_vertical_tangency_points
            (const Control_points& cp,
             const NT& t_min, const NT& t_max,
             std::list<Vertical_tangency_point>& vpts)
  {
    // \todo Handle the special case of degree two curves.

    // Check if we got to subdivision termination criteria.
    if (! can_refine (cp, t_min, t_max))
    {
      // It is not possible to further refine the approximation, so we stop
      // the recursive subdivision, and construct an output vertical tangency
      // point with can_refine = false. Note that we do not need a tight
      // bounding box here, since we cannot refine it anyway, so we take cp's
      // bbox (which certainly contains the intersection point).
      Bez_point_bound   bound (Bez_point_bound::VERTICAL_TANGENCY_PT,
                               cp, t_min, t_max,
                               false);
      Bez_point_bbox    vpt_bbox;

      construct_bbox (cp, vpt_bbox);
      vpts.push_back (Vertical_tangency_point (bound,
                                               vpt_bbox));

      return;
    }

    // If the control polygon is x-monotone, the curve does not contain any
    // vertical tangency points.
    if (_is_x_monotone (cp))
    {
      return;
    }

    // Check whether the control polygon is y-monotone.
    if (! _is_y_monotone (cp))
    {
      // Use de Casteljau's algorithm and subdivide the control polygon into
      // two, in order to obtain y-monotone subcurves.
      Control_points    cp_a, cp_b;
      const NT          t_mid = (t_min + t_max) / 2;

      bisect_control_polygon_2 (cp.begin(), cp.end(),
                                std::back_inserter (cp_a),
                                std::front_inserter (cp_b));

      // Check the case where t_mid is a vertical (rational) tangency
      // point, by checking whether the first vector of right is vertical.
      if (f_compare_x (cp_b[0], cp_b[1]) == EQUAL)
      {
        // The subdivision point at t_mid is a vertical tangency point with
        // rational coordinates.
        const typename Kernel::Point_2  &vpt = cp_b[0];

        Bez_point_bound   bound (Bez_point_bound::RATIONAL_PT,
                                 cp_b, t_mid, t_mid, true);
        Bez_point_bbox    bbox (vpt.x(), vpt.x(), vpt.y(), vpt.y());

        vpts.push_back (Vertical_tangency_point (bound,
                                                 bbox));
        return;
      }

      // Recursively compute the vertical tangency points of the two
      // subcurves.
      _compute_vertical_tangency_points (cp_a, t_min, t_mid,
                                         vpts);

      _compute_vertical_tangency_points (cp_b, t_mid, t_max,
                                         vpts);

      return;
    }

    // If we reached here, the control polygon of the curve is y-monotone,
    // but not x-monotone. Check whether it is convex. If it is, we know the
    // curve contains exactly one vertical tangency point.
    if (is_convex_2 (cp.begin(), cp.end(), m_kernel))
    {
      // We managed to isolate a single vertical tangency point.
      Bez_point_bound   bound (Bez_point_bound::VERTICAL_TANGENCY_PT,
                               cp, t_min, t_max,
                               true);    // It is possible to refine further.
      Bez_point_bbox    vpt_bbox;

      construct_bbox (cp, vpt_bbox);
      vpts.push_back (Vertical_tangency_point (bound,
                                               vpt_bbox));

      return;
    }

    // If we reached here, the curve contains more than a single vertical
    // tangency point. We therefore bisect it and continue recursively.
    Control_points    cp_a, cp_b;
    const NT          t_mid = (t_min + t_max) / 2;

    bisect_control_polygon_2 (cp.begin(), cp.end(),
                              std::back_inserter(cp_a),
                              std::front_inserter(cp_b));

    // Check the case where t_mid is a vertical (rational) tangency
    // point, by checking whether the first vector of right is vertical.
    if (f_compare_x (cp_b[0], cp_b[1]) == EQUAL)
    {
      // The subdivision point at t_mid is a vertical tangency point with
      // rational coordinates.
      const typename Kernel::Point_2  &vpt = cp_b[0];

      Bez_point_bound   bound (Bez_point_bound::RATIONAL_PT,
                               cp_b, t_mid, t_mid, true);
      Bez_point_bbox    bbox (vpt.x(), vpt.x(), vpt.y(), vpt.y());

      vpts.push_back (Vertical_tangency_point (bound,
                                               bbox));
      return;
    }

    // Recursively compute the vertical tangency points on the two subcurves.
    _compute_vertical_tangency_points (cp_a, t_min, t_mid,
                                       vpts);

    _compute_vertical_tangency_points (cp_b, t_mid, t_max,
                                       vpts);

    return;
  }
};

} //namespace CGAL

#endif //CGAL_BEZIER_BOUNDING_RATIONAL_TRAITS_H
