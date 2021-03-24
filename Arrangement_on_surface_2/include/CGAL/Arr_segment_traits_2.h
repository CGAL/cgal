// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Ron Wein          <wein@post.tau.ac.il>
//            Efi Fogel         <efif@gmail.com>
//            Waqar Khan        <wkhan@mpi-inf.mpg.de>
//            Simon Giraudot    <simon.giraudot@geometryfactory.com>

#ifndef CGAL_ARR_SEGMENT_TRAITS_2_H
#define CGAL_ARR_SEGMENT_TRAITS_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>

/*! \file
 * The segment traits-class for the arrangement package.
 */

#include <fstream>

#include <boost/variant.hpp>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/tags.h>
#include <CGAL/intersections.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_enums.h>
#include <CGAL/Arr_geometry_traits/Segment_assertions.h>

namespace CGAL {

template <typename Kernel_ = Exact_predicates_exact_constructions_kernel>
class Arr_segment_2;

/*! \class A traits class for maintaining an arrangement of segments, avoiding
 * cascading of computations as much as possible.
 *
 * The class is derived from the parameterized kernel to extend the traits
 * with all the types and operations supported by the kernel. This makes it
 * possible to use the traits class for data structures that extend the
 * Arrangement_2 type and require objects and operations supported by the
 * kernel, but not defined in this derived class.
 */
template <typename Kernel_ = Exact_predicates_exact_constructions_kernel>
class Arr_segment_traits_2 : public Kernel_ {
  friend class Arr_segment_2<Kernel_>;

public:
  typedef Kernel_                         Kernel;
  typedef typename Kernel::FT             FT;

  typedef typename Algebraic_structure_traits<FT>::Is_exact
                                          Has_exact_division;

  // Category tags:
  typedef Tag_true                        Has_left_category;
  typedef Tag_true                        Has_merge_category;
  typedef Tag_false                       Has_do_intersect_category;

  typedef Arr_oblivious_side_tag          Left_side_category;
  typedef Arr_oblivious_side_tag          Bottom_side_category;
  typedef Arr_oblivious_side_tag          Top_side_category;
  typedef Arr_oblivious_side_tag          Right_side_category;

  typedef typename Kernel::Line_2         Line_2;
  typedef CGAL::Segment_assertions<Arr_segment_traits_2<Kernel> >
                                          Segment_assertions;

  /*! \class Representation of a segment with cached data.
   */
  class _Segment_cached_2 {
  public:
    typedef typename Kernel::Line_2                Line_2;
    typedef typename Kernel::Segment_2             Segment_2;
    typedef typename Kernel::Point_2               Point_2;

  protected:
    mutable Line_2 m_l;         // the line that supports the segment.
    Point_2 m_ps;               // the source point of the segment.
    Point_2 m_pt;               // the target point of the segment.
    bool m_is_directed_right;   // is (lexicographically) directed left to right.
    mutable bool m_is_vert;     // is this a vertical segment.
    mutable bool m_is_computed; // is the support line computed.
    bool m_is_degen;            // is the segment degenerate (a single point).

  public:

    /// \name Creation
    //@{

    /*! Construct default. */
    _Segment_cached_2();

    /*! Construct a segment from a Kernel segment.
     * \param seg the segment.
     * \pre the segment is not degenerate.
     */
    _Segment_cached_2(const Segment_2& seg);

    /*! Construct a segment from two endpoints.
     * \param source the source point.
     * \param target the target point.
     * \param `source` and `target` are not equal.
     */
    _Segment_cached_2(const Point_2& source, const Point_2& target);

    /*! Construct a segment from two endpoints on a supporting line.
     * \param line the supporting line.
     * \param source the source point.
     * \param target the target point.
     * \pre `source` and `target` are not equal and both lie on `line`.
     */
    _Segment_cached_2(const Line_2& line,
                      const Point_2& source, const Point_2& target);

    /*! Construct a segment from all fields.
     * \param line the supporting line.
     * \param source the source point.
     * \param target the target point.
     * \param is_directed_right is (lexicographically) directed left to right.
     * \param is_vert is the segment vertical.
     * \param is_degen is the segment degenerate (a single point).
     */
    _Segment_cached_2(const Line_2& line,
                      const Point_2& source, const Point_2& target,
                      bool is_directed_right, bool is_vert, bool is_degen);

    /*! Assign.
     * \param seg the source segment to copy from
     * \pre the segment is not degenerate.
     */
    const _Segment_cached_2& operator=(const Segment_2& seg);

    //@}

    /// \name Accessors
    //@{

    /*! Obtain the supporting line.
     * \return the supporting line.
     */
    const Line_2& line() const;

    /*! Obtain the segment source.
     * \return the segment source.
     */
    const Point_2& source() const;

    /*! Obtain the segment target.
     * \return the segment target.
     */
    const Point_2& target() const;

    /*! Determine whether the curve is vertical.
     * \return a Boolean flag indicating whether the curve is vertical.
     */
    bool is_vertical() const;

    /*! Determine whether the curve is degenerate.
     * return a Boolean flag indicating whether the curve is degenerate.
     */
    bool is_degenerate() const;

    /*! Determine whether the curve is lexicographically directed from left to
     * right.
     * \return a Boolean flag indicating whether the curve is lexicographically
     *         directed from left to right.
     */
    bool is_directed_right() const;

    /*! Obtain the (lexicographically) left endpoint.
     * \return the (lexicographically) left endpoint.
     */
    const Point_2& left() const;

    /*! Obtain the (lexicographically) right endpoint.
     * \return the (lexicographically) right endpoint.
     */
    const Point_2& right() const;

    //@}

    /// \name Modifiers
    //@{

    /*! Set the (lexicographically) left endpoint.
     * \param p the point to set.
     * \pre p lies on the supporting line to the left of the right endpoint.
     */
    void set_left(const Point_2& p);

    /*! Set the (lexicographically) right endpoint.
     * \param p the point to set.
     * \pre p lies on the supporting line to the right of the left endpoint.
     */
    void set_right(const Point_2& p);

    //@}

    /// \name Deprecated
    //@{

    /*! Determine whether the given point is in the x-range of the segment.
     * \param p the query point.
     * \return (true) is in the x-range of the segment; (false) if it is not.
     */
    CGAL_DEPRECATED bool is_in_x_range(const Point_2& p) const;

    /*! Determine whether the given point is in the y-range of the segment.
     * \param p the query point.
     * \return (true) is in the y-range of the segment; (false) if it is not.
     */
    CGAL_DEPRECATED bool is_in_y_range(const Point_2& p) const;

    //@}
  };

public:
  // Traits objects
  typedef typename Kernel::Point_2        Point_2;
  typedef Arr_segment_2<Kernel>           X_monotone_curve_2;
  typedef Arr_segment_2<Kernel>           Curve_2;
  typedef unsigned int                    Multiplicity;

public:
  /*! Construct default. */
  Arr_segment_traits_2() {}

  /// \name Basic functor definitions.
  //@{

  class Compare_x_2 {
  protected:
    typedef Arr_segment_traits_2<Kernel>        Traits;

    //! The traits (in case it has state).
    const Traits& m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Compare_x_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_segment_traits_2<Kernel>;

  public:
    /*! Compare the x-coordinates of two points.
     * \param p1 the first point.
     * \param p2 the second point.
     * \return LARGER if x(p1) > x(p2);
     *         SMALLER if x(p1) < x(p2);
     *         EQUAL if x(p1) = x(p2).
     */
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const
    {
      const Kernel& kernel = m_traits;
      return (kernel.compare_x_2_object()(p1, p2));
    }
  };

  /*! Obtain a Compare_x_2 functor object. */
  Compare_x_2 compare_x_2_object() const { return Compare_x_2(*this); }

  class Compare_xy_2 {
  protected:
    typedef Arr_segment_traits_2<Kernel>        Traits;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Compare_xy_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_segment_traits_2<Kernel>;

  public:
    /*! Compare two points lexicographically: by x, then by y.
     * \param p1 the first point.
     * \param p2 the second point.
     * \return LARGER if x(p1) > x(p2), or if x(p1) = x(p2) and y(p1) > y(p2);
     *         SMALLER if x(p1) < x(p2), or if x(p1) = x(p2) and y(p1) < y(p2);
     *         EQUAL if the two points are equal.
     */
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const
    {
      const Kernel& kernel = m_traits;
      return (kernel.compare_xy_2_object()(p1, p2));
    }
  };

  /*! Obtain a Compare_xy_2 functor object. */
  Compare_xy_2 compare_xy_2_object() const { return Compare_xy_2(*this); }

  class Construct_min_vertex_2 {
  public:
    /*! Obtain the left endpoint of the x-monotone curve (segment).
     * \param cv the curve.
     * \return the left endpoint.
     */
    const Point_2& operator()(const X_monotone_curve_2& cv) const
    { return (cv.left()); }
  };

  /*! Obtain a Construct_min_vertex_2 functor object. */
  Construct_min_vertex_2 construct_min_vertex_2_object() const
  { return Construct_min_vertex_2(); }

  class Construct_max_vertex_2 {
  public:
    /*! Obtain the right endpoint of the x-monotone curve (segment).
     * \param cv the curve.
     * \return the right endpoint.
     */
    const Point_2& operator()(const X_monotone_curve_2& cv) const
    { return (cv.right()); }
  };

  /*! Obtain a Construct_max_vertex_2 functor object. */
  Construct_max_vertex_2 construct_max_vertex_2_object() const
  { return Construct_max_vertex_2(); }

  class Is_vertical_2 {
  public:
    /*! Check whether the given x-monotone curve is a vertical segment.
     * \param cv the curve.
     * \return (true) if the curve is a vertical segment; (false) otherwise.
     */
    bool operator()(const X_monotone_curve_2& cv) const
    { return (cv.is_vertical()); }
  };

  /*! Obtain an Is_vertical_2 functor object. */
  Is_vertical_2 is_vertical_2_object () const { return Is_vertical_2(); }

  class Compare_y_at_x_2 {
  protected:
    typedef Arr_segment_traits_2<Kernel>        Traits;

    /*! the traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Compare_y_at_x_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_segment_traits_2<Kernel>;

  public:
    /*! Return the location of the given point with respect to the input curve.
     * \param cv the curve.
     * \param p the point.
     * \pre p is in the x-range of cv.
     * \return SMALLER if y(p) < cv(x(p)), i.e. the point is below the curve;
     *         LARGER if y(p) > cv(x(p)), i.e. the point is above the curve;
     *         EQUAL if p lies on the curve.
     */
    Comparison_result operator()(const Point_2& p,
                                 const X_monotone_curve_2& cv) const
    {
      CGAL_precondition(m_traits.is_in_x_range_2_object()(cv, p));

      const Kernel& kernel = m_traits;

      if (! cv.is_vertical()) {
        // Compare p with the segment supporting line.
        CGAL_assertion_code(auto cmp_x = kernel.compare_x_2_object());
        CGAL_assertion(cmp_x(cv.left(), cv.right()) == SMALLER);
        return kernel.orientation_2_object()(cv.left(), cv.right(), p);
      }

      // Compare with the vertical segment endpoints.
      auto compare_y = kernel.compare_y_2_object();
      Comparison_result res1 = compare_y(p, cv.left());
      Comparison_result res2 = compare_y(p, cv.right());
      return (res1 == res2) ? res1 : EQUAL;
    }
  };

  /*! Obtain a Compare_y_at_x_2 functor object. */
  Compare_y_at_x_2 compare_y_at_x_2_object() const
  { return Compare_y_at_x_2(*this); }

  class Compare_y_at_x_left_2 {
  protected:
    typedef Arr_segment_traits_2<Kernel>        Traits;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Compare_y_at_x_left_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_segment_traits_2<Kernel>;

  public:
    /*! Compare the y value of two x-monotone curves immediately to the left
     * of their intersection point.
     * \param cv1 the first curve.
     * \param cv2 the second curve.
     * \param p the intersection point.
     * \pre the point p lies on both curves, and both of them must be also be
     *      defined (lexicographically) to its left.
     * \return the relative position of cv1 with respect to cv2 immediately to
     *         the left of p: SMALLER, LARGER or EQUAL.
     */
    Comparison_result operator()(const X_monotone_curve_2& cv1,
                                 const X_monotone_curve_2& cv2,
                                 const Point_2& CGAL_assertion_code(p)) const
    {
      const Kernel& kernel = m_traits;

      // Make sure that p lies on both curves, and that both are defined to its
      // left (so their left endpoint is lexicographically smaller than p).
      CGAL_precondition_code(auto compare_xy = kernel.compare_xy_2_object());

      CGAL_precondition((m_traits.compare_y_at_x_2_object()(p, cv1) == EQUAL) &&
                        (m_traits.compare_y_at_x_2_object()(p, cv2) == EQUAL));

      CGAL_precondition(compare_xy(cv1.left(), p) == SMALLER &&
                        compare_xy(cv2.left(), p) == SMALLER);

      // Compare the slopes of the two segments to determine their relative
      // position immediately to the left of q.
      // Notice we use the supporting lines in order to compare the slopes,
      // and that we swap the order of the curves in order to obtain the
      // correct result to the left of p.
      return (kernel.compare_slope_2_object()(cv2.line(), cv1.line()));
    }
  };

  /*! Obtain a Compare_y_at_x_left_2 functor object. */
  Compare_y_at_x_left_2 compare_y_at_x_left_2_object() const
  { return Compare_y_at_x_left_2(*this); }

  class Compare_y_at_x_right_2 {
  protected:
    typedef Arr_segment_traits_2<Kernel>        Traits;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Compare_y_at_x_right_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_segment_traits_2<Kernel>;

  public:
    /*! Compare the y value of two x-monotone curves immediately to the right
     * of their intersection point.
     * \param cv1 the first curve.
     * \param cv2 the second curve.
     * \param p the intersection point.
     * \pre the point p lies on both curves, and both of them must be also be
     *      defined (lexicographically) to its right.
     * \return the relative position of cv1 with respect to cv2 immediately to
     *         the right of p: SMALLER, LARGER or EQUAL.
     */
    Comparison_result operator()(const X_monotone_curve_2& cv1,
                                 const X_monotone_curve_2& cv2,
                                 const Point_2& CGAL_assertion_code(p)) const
    {
      const Kernel& kernel = m_traits;

      // Make sure that p lies on both curves, and that both are defined to its
      // right (so their right endpoint is lexicographically larger than p).
      CGAL_precondition_code(auto compare_xy = kernel.compare_xy_2_object());

      CGAL_precondition((m_traits.compare_y_at_x_2_object()(p, cv1) == EQUAL) &&
                        (m_traits.compare_y_at_x_2_object()(p, cv2) == EQUAL));

      CGAL_precondition(compare_xy(cv1.right(), p) == LARGER &&
                        compare_xy(cv2.right(), p) == LARGER);

      // Compare the slopes of the two segments to determine their relative
      // position immediately to the left of q.
      // Notice we use the supporting lines in order to compare the slopes.
      return (kernel.compare_slope_2_object()(cv1.line(), cv2.line()));
    }
  };

  /*! Obtain a Compare_y_at_x_right_2 functor object. */
  Compare_y_at_x_right_2 compare_y_at_x_right_2_object() const
  { return Compare_y_at_x_right_2(*this); }

  class Equal_2 {
  protected:
    typedef Arr_segment_traits_2<Kernel>        Traits;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Equal_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_segment_traits_2<Kernel>;

  public:
    /*! Check whether the two x-monotone curves are the same (have the same
     * graph).
     * \param cv1 the first curve.
     * \param cv2 the second curve.
     * \return (true) if the two curves are the same; (false) otherwise.
     */
    bool operator()(const X_monotone_curve_2& cv1,
                    const X_monotone_curve_2& cv2) const
    {
      const Kernel& kernel = m_traits;
      typename Kernel::Equal_2  equal = kernel.equal_2_object();

      return (equal(cv1.left(), cv2.left()) &&
              equal(cv1.right(), cv2.right()));
    }

    /*! Determine whether the two points are the same.
     * \param p1 the first point.
     * \param p2 the second point.
     * \return (true) if the two point are the same; (false) otherwise.
     */
    bool operator()(const Point_2& p1, const Point_2& p2) const
    {
      const Kernel& kernel = m_traits;
      return (kernel.equal_2_object()(p1, p2));
    }
  };

  /*! Obtain an Equal_2 functor object. */
  Equal_2 equal_2_object() const { return Equal_2(*this); }

  //@}

  //! \name Intersections, subdivisions, and mergings
  //@{

  /*! \class Make_x_monotone_2
   * A functor for subdividing a curve into x-monotone curves.
   */
  class Make_x_monotone_2 {
  public:
    /*! Subdivide a given curve into x-monotone subcurves and insert them into
     * a given output iterator. As segments are always x_monotone a single
     * object is inserted.
     * \param cv the curve.
     * \param oi the output iterator for the result. Its dereference type is a
     *           variant that wraps a \c Point_2 or an \c X_monotone_curve_2
     *           objects.
     * \return the past-the-end output iterator.
     */
    template <typename OutputIterator>
    OutputIterator operator()(const Curve_2& cv, OutputIterator oi) const
    {
      // Wrap the segment with a variant.
      typedef boost::variant<Point_2, X_monotone_curve_2>
        Make_x_monotone_result;
      *oi++ = Make_x_monotone_result(cv);
      return oi;
    }
  };

  /*! Obtain a Make_x_monotone_2 functor object. */
  Make_x_monotone_2 make_x_monotone_2_object() const
  { return Make_x_monotone_2(); }

  class Split_2 {
  protected:
    typedef Arr_segment_traits_2<Kernel>        Traits;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Split_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_segment_traits_2<Kernel>;

  public:
    /*! Split a given x-monotone curve at a given point into two sub-curves.
     * \param cv the curve to split
     * \param p the split point.
     * \param c1 Output: the left resulting subcurve (p is its right endpoint).
     * \param c2 Output: the right resulting subcurve (p is its left endpoint).
     * \pre p lies on cv but is not one of its endpoints.
     */
    void operator()(const X_monotone_curve_2& cv, const Point_2& p,
                    X_monotone_curve_2& c1, X_monotone_curve_2& c2) const
    {
      // Make sure that p lies on the interior of the curve.
      CGAL_precondition_code(const Kernel& kernel = m_traits;
                             auto compare_xy = kernel.compare_xy_2_object());

      CGAL_precondition((m_traits.compare_y_at_x_2_object()(p, cv) == EQUAL) &&
                        compare_xy(cv.left(), p) == SMALLER &&
                        compare_xy(cv.right(), p) == LARGER);

      // Perform the split.
      c1 = cv;
      c1.set_right(p);

      c2 = cv;
      c2.set_left(p);
    }
  };

  /*! Obtain a Split_2 functor object. */
  Split_2 split_2_object() const { return Split_2(*this); }

  class Intersect_2 {
  protected:
    typedef Arr_segment_traits_2<Kernel>        Traits;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Construct
     * \param traits the traits (in case it has state)
     */
    Intersect_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_segment_traits_2<Kernel>;

    // Specialized do_intersect with many tests skipped because at
    // this point, we already know which point is left / right for
    // both segments
    bool do_intersect(const Point_2& A1, const Point_2& A2,
                      const Point_2& B1, const Point_2& B2) const
    {
      const Kernel& kernel = m_traits;
      auto compare_xy = kernel.compare_xy_2_object();
      namespace interx = CGAL::Intersections::internal;

      switch(make_certain(compare_xy(A1,B1))) {
       case SMALLER:
        switch(make_certain(compare_xy(A2,B1))) {
         case SMALLER: return false;
         case EQUAL: return true;
         default: // LARGER
          switch(make_certain(compare_xy(A2,B2))) {
           case SMALLER:
            return interx::seg_seg_do_intersect_crossing(A1,A2,B1,B2, kernel);
           case EQUAL: return true;
           default: // LARGER
            return interx::seg_seg_do_intersect_contained(A1,A2,B1,B2, kernel);
          }
        }
       case EQUAL: return true;
       default: // LARGER
        switch(make_certain(compare_xy(B2,A1))) {
         case SMALLER: return false;
         case EQUAL: return true;
         default: // LARGER
          switch(make_certain(compare_xy(B2,A2))) {
           case SMALLER:
            return interx::seg_seg_do_intersect_crossing(B1,B2,A1,A2, kernel);
           case EQUAL: return true;
           default: // LARGER
            return interx::seg_seg_do_intersect_contained(B1,B2,A1,A2, kernel);
          }
        }
      }
      CGAL_error();     // never reached
      return false;
    }

    /*! Determine whether the bounding boxes of two segments overlap
     */
    bool do_bboxes_overlap(const X_monotone_curve_2& cv1,
                           const X_monotone_curve_2& cv2) const
    {
      const Kernel& kernel = m_traits;
      auto construct_bbox = kernel.construct_bbox_2_object();
      auto bbox1 = construct_bbox(cv1.source()) + construct_bbox(cv1.target());
      auto bbox2 = construct_bbox(cv2.source()) + construct_bbox(cv2.target());
      return CGAL::do_overlap(bbox1, bbox2);
    }

  public:
    /*! Find the intersections of the two given curves and insert them into the
     * given output iterator. As two segments may intersect only once, only a
     * single intersection will be contained in the iterator.
     * \param cv1 the first curve.
     * \param cv2 the second curve.
     * \param oi the output iterator.
     * \return the past-the-end iterator.
     */
    template <typename OutputIterator>
    OutputIterator operator()(const X_monotone_curve_2& cv1,
                              const X_monotone_curve_2& cv2,
                              OutputIterator oi) const
    {
      typedef std::pair<Point_2, Multiplicity>          Intersection_point;
      typedef boost::variant<Intersection_point, X_monotone_curve_2>
                                                        Intersection_result;

      // Early ending with Bbox overlapping test
      if (! do_bboxes_overlap(cv1, cv2)) return oi;

      // Early ending with specialized do_intersect
      const Kernel& kernel = m_traits;
      if (! do_intersect(cv1.left(), cv1.right(), cv2.left(), cv2.right()))
        return oi;

      // An intersection is guaranteed.

      // Intersect the two supporting lines.
      auto res = kernel.intersect_2_object()(cv1.line(), cv2.line());
      CGAL_assertion(bool(res));

      // Check if we have a single intersection point.
      const Point_2* ip = boost::get<Point_2>(&*res);
      if (ip != nullptr) {
        CGAL_assertion(cv1.is_vertical() ?
                       m_traits.is_in_y_range_2_object()(cv1, *ip) :
                       m_traits.is_in_x_range_2_object()(cv1, *ip));
        CGAL_assertion(cv2.is_vertical() ?
                       m_traits.is_in_y_range_2_object()(cv2, *ip) :
                       m_traits.is_in_x_range_2_object()(cv2, *ip));
        Intersection_point ip_mult(*ip, 1);
        *oi++ = Intersection_result(ip_mult);
        return oi;
      }

      // In this case, the two supporting lines overlap.
      // The overlapping segment is therefore [p_l,p_r], where p_l is the
      // rightmost of the two left endpoints and p_r is the leftmost of the
      // two right endpoints.
      auto compare_xy = kernel.compare_xy_2_object();
      const Point_2& p_l = (compare_xy(cv1.left(), cv2.left()) == SMALLER) ?
        cv2.left() : cv1.left();
      const Point_2& p_r = (compare_xy(cv1.right(), cv2.right()) == SMALLER) ?
        cv1.right() : cv2.right();

      // Examine the resulting segment.
      const Comparison_result cmp_res = compare_xy(p_l, p_r);
      if (cmp_res == EQUAL) {
        // The two segment have the same supporting line, but they just share
        // a common endpoint. Thus we have an intersection point, but we leave
        // the multiplicity of this point undefined.
        Intersection_point ip_mult(p_r, 0);
        *oi++ = Intersection_result(ip_mult);
        return oi;
      }

      CGAL_assertion(cmp_res == SMALLER);
      // We have discovered an overlapping segment:
      if (cv1.is_directed_right() == cv2.is_directed_right()) {
        // cv1 and cv2 have the same directions, maintain this direction
        // in the overlap segment
        if (cv1.is_directed_right()) {
          X_monotone_curve_2 overlap_seg(cv1.line(), p_l, p_r);
          *oi++ = Intersection_result(overlap_seg);
          return oi;
        }
        X_monotone_curve_2 overlap_seg(cv1.line(), p_r, p_l);
        *oi++ = Intersection_result(overlap_seg);
        return oi;
      }
      // cv1 and cv2 have opposite directions, the overlap segment
      // will be directed from left to right
      X_monotone_curve_2 overlap_seg(cv1.line(), p_l, p_r);
      *oi++ = Intersection_result(overlap_seg);
      return oi;
    }
  };

  /*! Obtain an Intersect_2 functor object. */
  Intersect_2 intersect_2_object() const { return Intersect_2(*this); }

  class Are_mergeable_2 {
  protected:
    typedef Arr_segment_traits_2<Kernel>        Traits;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Are_mergeable_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_segment_traits_2<Kernel>;

  public:
    /*! Check whether it is possible to merge two given x-monotone curves.
     * \param cv1 the first curve.
     * \param cv2 the second curve.
     * \return (true) if the two curves are mergeable, that is, if they are
     *         supported by the same line; (false) otherwise.
     * \pre cv1 and cv2 share a common endpoint.
     */
    bool operator()(const X_monotone_curve_2& cv1,
                    const X_monotone_curve_2& cv2) const
    {
      const Kernel& kernel = m_traits;
      typename Kernel::Equal_2 equal = kernel.equal_2_object();
      if (! equal(cv1.right(), cv2.left()) &&
          ! equal(cv2.right(), cv1.left()))
        return false;

      // Check whether the two curves have the same supporting line.
      return (equal(cv1.line(), cv2.line()) ||
              equal(cv1.line(),
                    kernel.construct_opposite_line_2_object()(cv2.line())));
    }
  };

  /*! Obtain an Are_mergeable_2 functor object. */
  Are_mergeable_2 are_mergeable_2_object() const
  { return Are_mergeable_2(*this); }

  /*! \class Merge_2
   * A functor that merges two x-monotone arcs into one.
   */
  class Merge_2 {
  protected:
    typedef Arr_segment_traits_2<Kernel>        Traits;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Merge_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_segment_traits_2<Kernel>;

  public:
    /*! Merge two given x-monotone curves into a single curve (segment).
     * \param cv1 the first curve.
     * \param cv2 the second curve.
     * \param c Output: the merged curve.
     * \pre the two curves are mergeable.
     */
    void operator()(const X_monotone_curve_2& cv1,
                    const X_monotone_curve_2& cv2,
                    X_monotone_curve_2& c) const
    {
      CGAL_precondition(m_traits.are_mergeable_2_object()(cv1, cv2));

      const Kernel& kernel = m_traits;
      auto equal = kernel.equal_2_object();

      // Check which curve extends to the right of the other.
      if (equal(cv1.right(), cv2.left())) {
        // cv2 extends cv1 to the right.
        c = cv1;
        c.set_right(cv2.right());
        return;
      }

      CGAL_precondition(equal(cv2.right(), cv1.left()));

      // cv1 extends cv2 to the right.
      c = cv2;
      c.set_right(cv1.right());
    }
  };

  /*! Obtain a Merge_2 functor object. */
  Merge_2 merge_2_object() const { return Merge_2(*this); }
  //@}

  /// \name Functor definitions for the landmarks point-location strategy.
  //@{
  typedef double                          Approximate_number_type;

  class Approximate_2 {
  public:
    /*! Obtain an approximation of a point coordinate.
     * \param p the exact point.
     * \param i the coordinate index (either 0 or 1).
     * \pre i is either 0 or 1.
     * \return An approximation of p's x-coordinate (if i == 0), or an
     *         approximation of p's y-coordinate (if i == 1).
     */
    Approximate_number_type operator()(const Point_2& p, int i) const
    {
      CGAL_precondition((i == 0) || (i == 1));
      return (i == 0) ? (CGAL::to_double(p.x())) : (CGAL::to_double(p.y()));
    }
  };

  /*! Obtain an Approximate_2 functor object. */
  Approximate_2 approximate_2_object() const { return Approximate_2(); }

  class Construct_x_monotone_curve_2 {
  protected:
    typedef Arr_segment_traits_2<Kernel>        Traits;

    //! The traits (in case it has state).
    const Traits& m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Construct_x_monotone_curve_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_segment_traits_2<Kernel>;

  public:
    typedef typename Kernel::Segment_2          Segment_2;

    /*! Obtain an x-monotone curve connecting two given endpoints.
     * \param source the first point.
     * \param target the second point.
     * \pre `source` and `target` must not be equal.
     * \return a segment connecting `source` and `target`.
     */
    X_monotone_curve_2 operator()(const Point_2& source,
                                  const Point_2& target) const
    {
      const Kernel& kernel = m_traits;
      auto line = kernel.construct_line_2_object()(source, target);
      Comparison_result res = kernel.compare_xy_2_object()(source, target);
      auto is_degen = (res == EQUAL);
      auto is_directed_right = (res == SMALLER);
      CGAL_precondition_msg(! is_degen,
                            "Cannot construct a degenerate segment.");
      auto is_vert = kernel.is_vertical_2_object()(line);
      return X_monotone_curve_2(line, source, target,
                                is_directed_right, is_vert, is_degen);
    }

    /*! Obtain an \f$x\f$-monotone curve given a Kernel segment.
     * \param seg the segment.
     * \return the \f$x\f$-monotone curve.
     * \pre the segment is not degenerate.
     * \return a segment that is the same as `seg`..
     */
    X_monotone_curve_2 operator()(const Segment_2& seg) const
    {
      const Kernel& kernel = m_traits;
      auto line = kernel.construct_line_2_object()(seg);
      auto vertex_ctr = kernel.construct_vertex_2_object();
      auto source = vertex_ctr(seg, 0);
      auto target = vertex_ctr(seg, 1);
      Comparison_result res = kernel.compare_xy_2_object()(source, target);
      auto is_degen = (res == EQUAL);
      auto is_directed_right = (res == SMALLER);
      CGAL_precondition_msg(! is_degen,
                            "Cannot construct a degenerate segment.");
      auto is_vert = kernel.is_vertical_2_object()(seg);
      return X_monotone_curve_2(line, source, target,
                                is_directed_right, is_vert, is_degen);
    }

    /*! Obtain an \f$x\f$-monotone curve given two endpoints and the supporting
     * line.
     * \param line the supporting line.
     * \param  the source point.
     * \param  the target point.
     * \pre `ps` and `pt` are not equal and both lie on `line`.
     */
    X_monotone_curve_2 operator()(const Line_2& line,
                                  const Point_2& source,
                                  const Point_2& target) const
    {
      const Kernel& kernel = m_traits;
      CGAL_precondition
        (Segment_assertions::_assert_is_point_on(source, line,
                                                 Has_exact_division()) &&
         Segment_assertions::_assert_is_point_on(target, line,
                                                 Has_exact_division()));
      auto is_vert = kernel.is_vertical_2_object()(line);
      Comparison_result res = kernel.compare_xy_2_object()(source, target);
      auto is_degen = (res == EQUAL);
      auto is_directed_right = (res == SMALLER);
      CGAL_precondition_msg(! is_degen,
                            "Cannot construct a degenerate segment.");
      return X_monotone_curve_2(line, source, target,
                                is_directed_right, is_vert, is_degen);
    }
  };

  /*! Obtain a Construct_x_monotone_curve_2 functor object. */
  Construct_x_monotone_curve_2 construct_x_monotone_curve_2_object() const
  { return Construct_x_monotone_curve_2(*this); }
  //@}

  /// \name Functor definitions for the Boolean set-operation traits.
  //@{

  class Trim_2 {
   protected:
    typedef Arr_segment_traits_2<Kernel>        Traits;

    /*! The traits (in case it has state). */
    const Traits& m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Trim_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_segment_traits_2<Kernel>;

    /*! Obtain a trimmed version of a line.
     * \param xseg the x-monotone segment.
     * \param src the new start endpoint.
     * \param tgt the new end endpoint.
     * \return the trimmed x-monotone segment.
     * \pre src != tgt
     * \pre both points must lie on segment
     */
  public:
    X_monotone_curve_2 operator()(const X_monotone_curve_2& xcv,
                                  const Point_2& src,
                                  const Point_2& tgt)const
    {
      CGAL_precondition_code(Equal_2 equal = m_traits.equal_2_object());
      CGAL_precondition_code(Compare_y_at_x_2 compare_y_at_x =
                             m_traits.compare_y_at_x_2_object());
      Compare_x_2 compare_x_2 = m_traits.compare_x_2_object();

      // check whether source and taget are two distinct points and they lie
      // on the line.
      CGAL_precondition(!equal(src, tgt));
      CGAL_precondition(compare_y_at_x(src, xcv) == EQUAL);
      CGAL_precondition(compare_y_at_x(tgt, xcv) == EQUAL);

      // exchange src and tgt IF they do not conform with the direction
      X_monotone_curve_2 trimmed_segment;
      if (xcv.is_directed_right() && compare_x_2(src, tgt) == LARGER)
        trimmed_segment = X_monotone_curve_2(tgt, src);
      else if (! xcv.is_directed_right() && (compare_x_2(src, tgt) == SMALLER))
        trimmed_segment = X_monotone_curve_2(tgt, src);
      else trimmed_segment = X_monotone_curve_2(src, tgt);
      return trimmed_segment;
    }
  };

  /*! Obtain a Trim_2 functor object */
  Trim_2 trim_2_object() const { return Trim_2(*this); }

  class Compare_endpoints_xy_2 {
  public:
    /*! Compare the endpoints of an $x$-monotone curve lexicographically.
     * (assuming the curve has a designated source and target points).
     * \param cv the curve.
     * \return SMALLER if the curve is directed right;
     *         LARGER if the curve is directed left.
     */
    Comparison_result operator()(const X_monotone_curve_2& cv) const
    { return (cv.is_directed_right()) ? (SMALLER) : (LARGER); }
  };

  /*! Obtain a Compare_endpoints_xy_2 functor object. */
  Compare_endpoints_xy_2 compare_endpoints_xy_2_object() const
  { return Compare_endpoints_xy_2(); }

  class Construct_opposite_2 {
  public:
    /*! Construct an opposite x-monotone (with swapped source and target).
     * \param cv the curve.
     * \return the opposite curve.
     */
    X_monotone_curve_2 operator()(const X_monotone_curve_2& cv) const
    { return (cv.flip()); }
  };

  /*! Obtain a Construct_opposite_2 functor object. */
  Construct_opposite_2 construct_opposite_2_object() const
  { return Construct_opposite_2(); }
  //@}

  /// \name Utilities.
  //@{

  class Is_in_x_range_2 {
  protected:
    typedef Arr_segment_traits_2<Kernel>        Traits;

    //! The traits (in case it has state).
    const Traits& m_traits;

    /*! Construct
     * \param traits the traits (in case it has state)
     */
    Is_in_x_range_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_segment_traits_2<Kernel>;

  public:
    /*! Determine whether a given point is in the \f$x\f$-range of a given
     * segment.
     * \param cv the segment.
     * \param p the point.
     * \return true if p is in the \f$x\f$-range of cv; false otherwise.
     */
    bool operator()(const X_monotone_curve_2& cv, const Point_2& p) const
    {
      const Kernel& kernel = m_traits;
      auto compare_x = kernel.compare_x_2_object();
      Comparison_result res1 = compare_x(p, cv.left());

      if (res1 == SMALLER) return false;
      else if (res1 == EQUAL) return true;

      Comparison_result res2 = compare_x(p, cv.right());
      return (res2 != LARGER);
    }
  };

  /*! Obtain an Is_in_x_range_2 functor object */
  Is_in_x_range_2 is_in_x_range_2_object() const
  { return Is_in_x_range_2(*this); }

  class Is_in_y_range_2 {
  protected:
    typedef Arr_segment_traits_2<Kernel>        Traits;

    //! The traits (in case it has state).
    const Traits& m_traits;

    /*! Construct
     * \param traits the traits (in case it has state)
     */
    Is_in_y_range_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_segment_traits_2<Kernel>;

  public:
    /*! Determine whether a given point is in the \f$y\f$-range of a given
     * segment.
     * \param cv the segment.
     * \param p the point.
     * \return true if p is in the \f$y\f$-range of cv; false otherwise.
     */
    bool operator()(const X_monotone_curve_2& cv, const Point_2& p) const
    {
      const Kernel& kernel = m_traits;
      auto compare_y = kernel.compare_y_2_object();
      Comparison_result res1 = compare_y(p, cv.left());

      if (res1 == SMALLER) return false;
      else if (res1 == EQUAL) return true;

      Comparison_result res2 = compare_y(p, cv.right());
      return (res2 != LARGER);
    }
  };

  /*! Obtain an Is_in_y_range_2 functor object */
  Is_in_y_range_2 is_in_y_range_2_object() const
  { return Is_in_y_range_2(*this); }

  //@}
};

// Creation

//! \brief constructs default.
template <typename Kernel>
Arr_segment_traits_2<Kernel>::_Segment_cached_2::_Segment_cached_2() :
  m_is_directed_right(false),
  m_is_vert(false),
  m_is_computed(false),
  m_is_degen(true)
{}

//! \brief constructs a segment from a Kernel segment.
template <typename Kernel>
Arr_segment_traits_2<Kernel>::
_Segment_cached_2::_Segment_cached_2(const Segment_2& seg) :
  m_is_vert(false),
  m_is_computed(false)
{
  Kernel kernel;
  auto vertex_ctr = kernel.construct_vertex_2_object();

  m_ps = vertex_ctr(seg, 0);
  m_pt = vertex_ctr(seg, 1);

  Comparison_result res = kernel.compare_xy_2_object()(m_ps, m_pt);
  m_is_degen = (res == EQUAL);
  m_is_directed_right = (res == SMALLER);

  CGAL_precondition_msg(! m_is_degen, "Cannot construct a degenerate segment.");
}

//! \brief Constructs a segment from two endpoints.
template <typename Kernel>
Arr_segment_traits_2<Kernel>::
_Segment_cached_2::_Segment_cached_2(const Point_2& source,
                                     const Point_2& target) :
  m_ps(source),
  m_pt(target),
  m_is_vert(false),
  m_is_computed(false)
{
  Kernel kernel;

  Comparison_result res = kernel.compare_xy_2_object()(m_ps, m_pt);
  m_is_degen = (res == EQUAL);
  m_is_directed_right = (res == SMALLER);

  CGAL_precondition_msg(! m_is_degen, "Cannot construct a degenerate segment.");
}

//! \brief constructs a segment from two endpoints on a supporting line.
template <typename Kernel>
Arr_segment_traits_2<Kernel>::
_Segment_cached_2::_Segment_cached_2(const Line_2& line,
                                     const Point_2& source,
                                     const Point_2& target) :
  m_l(line),
  m_ps(source),
  m_pt(target)
{
  Kernel kernel;

  CGAL_precondition
    (Segment_assertions::_assert_is_point_on(source, m_l,
                                             Has_exact_division()) &&
     Segment_assertions::_assert_is_point_on(target, m_l,
                                             Has_exact_division()));

  m_is_vert = kernel.is_vertical_2_object()(m_l);
  m_is_computed = true;

  Comparison_result res = kernel.compare_xy_2_object()(m_ps, m_pt);
  m_is_degen = (res == EQUAL);
  m_is_directed_right = (res == SMALLER);

  CGAL_precondition_msg(! m_is_degen, "Cannot construct a degenerate segment.");
}

//! \brief constructs a segment from all fields.
template <typename Kernel>
Arr_segment_traits_2<Kernel>::_Segment_cached_2::
_Segment_cached_2(const Line_2& line,
                  const Point_2& source, const Point_2& target,
                  bool is_directed_right, bool is_vert, bool is_degen) :
  m_l(line),
  m_ps(source),
  m_pt(target),
  m_is_directed_right(is_directed_right),
  m_is_vert(is_vert),
  m_is_computed(true),
  m_is_degen(is_degen)
{}

//! \brief assigns.
template <typename Kernel>
const typename Arr_segment_traits_2<Kernel>::_Segment_cached_2&
Arr_segment_traits_2<Kernel>::_Segment_cached_2::operator=(const Segment_2& seg)
{
  Kernel kernel;
  auto vertex_ctr = kernel.construct_vertex_2_object();

  m_ps = vertex_ctr(seg, 0);
  m_pt = vertex_ctr(seg, 1);

  Comparison_result res = kernel.compare_xy_2_object()(m_ps, m_pt);
  m_is_degen = (res == EQUAL);
  m_is_directed_right = (res == SMALLER);

  CGAL_precondition_msg(! m_is_degen, "Cannot construct a degenerate segment.");

  m_l = kernel.construct_line_2_object()(seg);
  m_is_vert = kernel.is_vertical_2_object()(seg);
  m_is_computed = true;

  return (*this);
}

// Accessors

//! \brief obtains the supporting line.
template <typename Kernel>
const typename Kernel::Line_2&
Arr_segment_traits_2<Kernel>::_Segment_cached_2::line() const
{
  if (!m_is_computed) {
    Kernel kernel;
    m_l = kernel.construct_line_2_object()(m_ps, m_pt);
    m_is_vert = kernel.is_vertical_2_object()(m_l);
    m_is_computed = true;
  }
  return m_l;
}

//! \brief determines whether the curve is vertical.
template <typename Kernel>
bool Arr_segment_traits_2<Kernel>::_Segment_cached_2::is_vertical() const
{
  // Force computation of line is orientation is still unknown
  if (! m_is_computed) line();
  CGAL_precondition(!m_is_degen);
  return m_is_vert;
}

//! \brief determines whether the curve is degenerate.
template <typename Kernel>
bool Arr_segment_traits_2<Kernel>::_Segment_cached_2::is_degenerate() const
{ return m_is_degen; }

/*! \brief determines whether the curve is lexicographically directed from
 * left to right
 */
template <typename Kernel>
bool Arr_segment_traits_2<Kernel>::_Segment_cached_2::is_directed_right() const
{ return m_is_directed_right; }

//! \brief obtain the segment source.
template <typename Kernel>
const typename Kernel::Point_2&
Arr_segment_traits_2<Kernel>::_Segment_cached_2::source() const { return m_ps; }

//! \brief obtains the segment target.
template <typename Kernel>
const typename Kernel::Point_2&
Arr_segment_traits_2<Kernel>::_Segment_cached_2::target() const { return m_pt; }

//! \brief obtains the (lexicographically) left endpoint.
template <typename Kernel>
const typename Kernel::Point_2&
Arr_segment_traits_2<Kernel>::_Segment_cached_2::left() const
{ return (m_is_directed_right ? m_ps : m_pt); }

//! \brief obtains the (lexicographically) right endpoint.
template <typename Kernel>
const typename Kernel::Point_2&
Arr_segment_traits_2<Kernel>::_Segment_cached_2::right() const
{ return (m_is_directed_right ? m_pt : m_ps); }

// Modifiers

//! \brief sets the (lexicographically) left endpoint.
template <typename Kernel>
void Arr_segment_traits_2<Kernel>::_Segment_cached_2::set_left(const Point_2& p)
{
  CGAL_precondition(! m_is_degen);
  CGAL_precondition_code(Kernel kernel);
  CGAL_precondition
    (Segment_assertions::_assert_is_point_on(p, m_l, Has_exact_division()) &&
     (kernel.compare_xy_2_object()(p, right()) == SMALLER));

  if (m_is_directed_right) m_ps = p;
  else m_pt = p;
}

//! \brief sets the (lexicographically) right endpoint.
template <typename Kernel>
void Arr_segment_traits_2<Kernel>::_Segment_cached_2::set_right(const Point_2& p)
{
  CGAL_precondition(! m_is_degen);
  CGAL_precondition_code(Kernel kernel);
  CGAL_precondition
    (Segment_assertions::_assert_is_point_on(p, m_l, Has_exact_division()) &&
     (kernel.compare_xy_2_object()(p, left()) == LARGER));

  if (m_is_directed_right) m_pt = p;
  else m_ps = p;
}

//! \brief determines whether the given point is in the x-range of the segment.
template <typename Kernel>
bool Arr_segment_traits_2<Kernel>::_Segment_cached_2::
is_in_x_range(const Point_2& p) const
{
  Kernel kernel;
  typename Kernel::Compare_x_2 compare_x = kernel.compare_x_2_object();
  const Comparison_result res1 = compare_x(p, left());

  if (res1 == SMALLER) return false;
  else if (res1 == EQUAL) return true;

  const Comparison_result res2 = compare_x(p, right());
  return (res2 != LARGER);
}

//! \brief determines whether the given point is in the y-range of the segment.
template <typename Kernel>
bool Arr_segment_traits_2<Kernel>::_Segment_cached_2::
is_in_y_range(const Point_2& p) const
{
  Kernel kernel;
  typename Kernel::Compare_y_2 compare_y = kernel.compare_y_2_object();
  const Comparison_result res1 = compare_y(p, left());

  if (res1 == SMALLER) return false;
  else if (res1 == EQUAL) return true;

  const Comparison_result res2 = compare_y(p, right());
  return (res2 != LARGER);
}

/*! \class A representation of a segment, as used by the Arr_segment_traits_2
 * traits-class.
 */
template <typename Kernel_>
class Arr_segment_2 : public Arr_segment_traits_2<Kernel_>::_Segment_cached_2 {
  typedef Kernel_                                                  Kernel;

  typedef typename Arr_segment_traits_2<Kernel>::_Segment_cached_2 Base;
  typedef typename Kernel::Segment_2                               Segment_2;
  typedef typename Kernel::Point_2                                 Point_2;
  typedef typename Kernel::Line_2                                  Line_2;

public:
  /*! Construct default. */
  Arr_segment_2();

  /*! Construct a segment from a "kernel" segment.
   * \param seg the segment.
   * \pre the segment is not degenerate.
   */
  Arr_segment_2(const Segment_2& seg);

  /*! Construct a segment from two endpoints.
   * \param source the source point.
   * \param target the target point.
   * \pre `source` and `target` are not equal.
   */
  Arr_segment_2(const Point_2& source, const Point_2& target);

  /*! Construct a segment from a line and two endpoints.
   * \param line the supporting line.
   * \param source the source point.
   * \param target the target point.
   * \pre Both `source` and `target` must be on `line`.
   * \pre `source` and `target` are not equal.
   */
  Arr_segment_2(const Line_2& line,
                const Point_2& source, const Point_2& target);

  /*! Construct a segment from all fields.
   * \param line the supporting line.
   * \param source the source point.
   * \param target the target point.
   * \param is_directed_right is (lexicographically) directed left to right.
   * \param is_vert is the segment vertical.
   * \param is_degen is the segment degenerate (a single point).
   */
  Arr_segment_2(const Line_2& line,
                const Point_2& source, const Point_2& target,
                bool is_directed_right, bool is_vert, bool is_degen);

  /*! Cast to a segment.
   */
  operator Segment_2() const;

  /*! Flip the segment (swap its source and target).
   */
  Arr_segment_2 flip() const;

  /*! Create a bounding box for the segment.
   */
  Bbox_2 bbox() const;
};

//! \brief constructs default.
template <typename Kernel>
Arr_segment_2<Kernel>::Arr_segment_2() : Base() {}

//! \brief constructs from a "kernel" segment.
template <typename Kernel>
Arr_segment_2<Kernel>::Arr_segment_2(const Segment_2& seg) : Base(seg) {}

//! \brief constructs a segment from two end-points.
template <typename Kernel>
Arr_segment_2<Kernel>::Arr_segment_2(const Point_2& source,
                                     const Point_2& target) :
  Base(source, target)
{}

//! \brief constructs a segment from a line and two end-points.
template <typename Kernel>
Arr_segment_2<Kernel>::Arr_segment_2(const Line_2& line,
                                     const Point_2& source,
                                     const Point_2& target) :
  Base(line,source, target)
{}

//! \brief constructs a segment from all fields.
template <typename Kernel>
Arr_segment_2<Kernel>::Arr_segment_2(const Line_2& line,
                                     const Point_2& source,
                                     const Point_2& target,
                                     bool is_directed_right,
                                     bool is_vert, bool is_degen) :
  Base(line, source, target, is_directed_right, is_vert, is_degen)
{}

//! \brief casts to a segment.
template <typename Kernel>
Arr_segment_2<Kernel>::operator typename Kernel::Segment_2() const
{
  Kernel kernel;
  auto seg_ctr = kernel.construct_segment_2_object();
  return seg_ctr(this->source(), this->target());
}

//! \brief flips the segment (swap its source and target).
template <typename Kernel>
Arr_segment_2<Kernel> Arr_segment_2<Kernel>::flip() const
{
  return Arr_segment_2(this->line(), this->target(), this->source(),
                       ! (this->is_directed_right()), this->is_vertical(),
                       this->is_degenerate());
}

//! \brief creates a bounding box for the segment.
template <typename Kernel>
Bbox_2 Arr_segment_2<Kernel>::bbox() const
{
  Kernel kernel;
  auto construct_bbox = kernel.construct_bbox_2_object();
  return construct_bbox(this->m_ps) + construct_bbox(this->m_pt);
}

/*! Exporter for the segment class used by the traits-class.
 */
template <typename Kernel, typename OutputStream>
OutputStream& operator<<(OutputStream& os, const Arr_segment_2<Kernel>& seg)
{
  os << static_cast<typename Kernel::Segment_2>(seg);
  return (os);
}

/*! Importer for the segment class used by the traits-class.
 */
template <typename Kernel, typename InputStream>
InputStream& operator>>(InputStream& is, Arr_segment_2<Kernel>& seg)
{
  typename Kernel::Segment_2   kernel_seg;
  is >> kernel_seg;
  seg = kernel_seg;
  return is;
}

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
