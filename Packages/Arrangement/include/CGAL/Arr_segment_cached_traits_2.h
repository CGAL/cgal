// Copyright (c) 1997  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Ron Wein          <wein@post.tau.ac.il>
//                 Efi Fogel         <efif@post.tau.ac.il>

#ifndef CGAL_ARR_SEGMENT_CACHED_TRAITS_2_H
#define CGAL_ARR_SEGMENT_CACHED_TRAITS_2_H

#include <CGAL/tags.h>
#include <CGAL/intersections.h>

#include <list>
#include <fstream>

CGAL_BEGIN_NAMESPACE

template <class Kernel_> class Segment_cached_2;

/*!
 * A traits class for maintaining an arrangement of segments, aoviding
 * cascading of computations as much as possible.
 */
template <class Kernel_>
class Arr_segment_cached_traits_2 : public Kernel_
{
public:
  typedef Kernel_                         Kernel;

  typedef Tag_true                        Has_left_category;
  typedef Tag_false                       Has_reflect_category;

  /*!
   * Representation of a segement with cached data.
   */
  class My_segment_cached_2
  {
    typedef typename Kernel_::Line_2                Line_2;
    typedef typename Kernel_::Segment_2             Segment_2;
    typedef typename Kernel_::Point_2               Point_2;

  protected:

    Line_2    line;             // The line that supports the segment.
    Point_2   ps, pt;           // The source a target points.
    bool      is_vert;          // Is this a vertical segment.

    /*!
     * Default constructor.
     */
    My_segment_cached_2 () :
      is_vert(false)
    {}

    /*!
     * Constructor from a segment.
     * \param seg The segment.
     */
    My_segment_cached_2 (const Segment_2 & seg)
    {
      Kernel_   kernel;
      
      line = kernel.construct_line_2_object()(seg);
      is_vert = kernel.is_vertical_2_object()(seg);
      
      typename Kernel_::Construct_vertex_2 
        construct_vertex = kernel.construct_vertex_2_object();

      ps = construct_vertex(seg, 0);
      pt = construct_vertex(seg, 1);
    }

    /*!
     * Construct a segment from two end-points.
     * \param source The source point.
     * \param target The target point.
     */
    My_segment_cached_2 (const Point_2 & source, const Point_2 & target) :
      ps (source),
      pt (target)
    {
      Kernel_   kernel;
      
      line = kernel.construct_line_2_object()(source, target);
      is_vert = kernel.is_vertical_2_object()(line);
    }

    /*!
     * Construct a segment from two end-points on a supporting line.
     * \param l The supporting line.
     * \param source The source point.
     * \param target The target point.
     */
    My_segment_cached_2 (const Line_2& l,
                         const Point_2 & source, const Point_2 & target) :
      line (l),
      ps (source),
      pt (target)
    {
      Kernel_   kernel;

      CGAL_precondition (kernel.has_on_2_object() (line, source) &&
                         kernel.has_on_2_object() (line, target));

      is_vert = kernel.is_vertical_2_object()(line);
    }

    /*!
     * Assignment operator.
     * \param seg the source segment to copy from
     */
    const My_segment_cached_2& operator= (const Segment_2 & seg)
    {
      Kernel_   kernel;
      
      line = kernel.construct_line_2_object()(seg);
      is_vert = kernel.is_vertical_2_object()(seg);
      
      typename Kernel_::Construct_vertex_2 
        construct_vertex = kernel.construct_vertex_2_object();

      ps = construct_vertex(seg, 0);
      pt = construct_vertex(seg, 1);
      return *this;                   
    }
    
    friend class Arr_segment_cached_traits_2;
  };

  // Traits objects
  typedef typename Kernel::Point_2        Point_2;
  typedef Segment_cached_2<Kernel>        X_monotone_curve_2;
  typedef Segment_cached_2<Kernel>        Curve_2;

  // Backward compatability    
  typedef Point_2                         Point;
  typedef X_monotone_curve_2              X_curve;
  typedef X_monotone_curve_2              Curve;

protected:

  // Functors:
  typedef typename Kernel::Less_x_2             Less_x_2;
  typedef typename Kernel::Equal_2              Equal_2;
  typedef typename Kernel::Compare_x_2          Compare_x_2;
  typedef typename Kernel::Compare_y_2          Compare_y_2;
  typedef typename Kernel::Compare_xy_2         Compare_xy_2;
  typedef typename Kernel::Compare_slope_2      Compare_slope_2;
  typedef typename Kernel::Has_on_2             Has_on_2;
    
  public:

#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_3
  using Kernel::compare_xy_2_object;
  using Kernel::compare_x_2_object;
  using Kernel::compare_y_2_object;
  using Kernel::compare_y_at_x_2_object;
  using Kernel::is_vertical_2_object;
  using Kernel::has_on_2_object;
  using Kernel::compare_slope_2_object;
  using Kernel::equal_2_object;
  using Kernel::intersect_2_object;
  using Kernel::construct_line_2_object;
#endif


  /*!
   * Default constructor.
   */
  Arr_segment_cached_traits_2() {}

  // Operations
  // ----------
    
  /*!
   * Compare the x-coordinates of two given points.
   * \param p1 The first point.
   * \param p2 The second point.
   * \return LARGER if x(p1) > x(p2); SMALLER if x(p1) < x(p2); or else EQUAL.
   */
  Comparison_result compare_x(const Point_2 & p1, const Point_2 & p2) const
  {
    return (compare_x_2_object()(p1, p2));
  }

  /*! 
   * Compares lexigoraphically the two points: by x, then by y.
   * \param p1 Te first point.
   * \param p2 The second point.
   * \return LARGER if x(p1) > x(p2), or if x(p1) = x(p2) and y(p1) > y(p2); 
   *         SMALLER if x(p1) < x(p2), or if x(p1) = x(p2) and y(p1) < y(p2);
   *         or else EQUAL.
   */
  Comparison_result compare_xy(const Point_2 & p1, const Point_2 & p2) const
  {
    return (compare_xy_2_object()(p1, p2));
  }

  /*!
   * Check whether the given curve is a vertical segment.
   * \param cv The curve.
   * \return (true) if the curve is vertical.
   */
  bool curve_is_vertical(const X_monotone_curve_2 & cv) const 
  {
    return (cv.is_vert);
  } 

  /*!
   * Check whether the given point is in the x-range of the given curve.
   * In out case, the curve is a segment [s, t], check whether x(s)<=x(q)<=x(t)
   * or whether x(t)<=x(q)<=x(s).
   * \param cv The curve.
   * \param q The point.
   * \return (true) if q is in the x-range of cv.
   */
  bool point_in_x_range(const X_monotone_curve_2 & cv, const Point_2 & q) const
  {
    Compare_x_2       compare_x = compare_x_2_object();
    Comparison_result res1 = compare_x (q, cv.ps);

    if (cv.is_vert) // Special check for vertical segments.
      return (res1 == EQUAL);

    Comparison_result res2 = compare_x (q, cv.pt);

    // We check if x(p) equals the x value of one of the end-points.
    // If not, we check whether one end-point is to p's left and the other is
    // to its right.
    return ((res1 == EQUAL) || (res2 == EQUAL) ||
            (res1 != res2));
  }

  /*!
   * Get the relative status of two curves at the x-coordinate of a given 
   * point.
   * \param cv1 The first curve.
   * \param cv2 The second curve.
   * \param q The point.
   * \pre The point q is in the x-range of the two curves.
   * \return LARGER if cv1(x(q)) > cv2(x(q)); SMALLER if cv1(x(q)) < cv2(x(q));
   *  or else EQUAL.
   */
  Comparison_result curves_compare_y_at_x(const X_monotone_curve_2 & cv1, 
                                          const X_monotone_curve_2 & cv2, 
                                          const Point_2 & q) const
  {
    CGAL_precondition(point_in_x_range(cv1, q));
    CGAL_precondition(point_in_x_range(cv2, q));

    // Special treatment for vertical segments:
    // In case one curve is a vertical segment, return EQUAL if it intersects
    // with the other segment, otherwise return LARGER or SMALLER.
    if (cv1.is_vert)
    {
      if (cv2.is_vert)
      {
        // Compare two vertical segments.
        Compare_y_2       compare_y = compare_y_2_object();
        Comparison_result res1 = compare_y (cv1.ps, cv1.pt);
        const Point_2 &   lower1 = (res1 == SMALLER) ? cv1.ps : cv1.pt;
        const Point_2 &   upper1 = (res1 == SMALLER) ? cv1.pt : cv1.ps;
        Comparison_result res2 = compare_y (cv2.ps, cv2.pt);
        const Point_2 &   lower2 = (res2 == SMALLER) ? cv2.ps : cv2.pt;
        const Point_2 &   upper2 = (res2 == SMALLER) ? cv2.pt : cv2.ps;

        if (compare_y(upper1, lower2) == SMALLER)
          // cv1 is entirely below cv2:
          return (SMALLER); 
        else if (compare_y(lower1, upper2) == LARGER)
          // cv1 is entirely above cv2:
          return (LARGER);
        else
          // cv1 intersects cv2:
          return (EQUAL);
      }

      // Only cv1 is vertical:
      Comparison_result res1 = compare_y_at_x_2_object()(cv1.ps,cv2.line);
      Comparison_result res2 = compare_y_at_x_2_object()(cv1.pt,cv2.line);

      if (res1 == res2)
      {
        CGAL_assertion(res1 != EQUAL);
        return (res1);
      }
      else
        return (EQUAL);
    }
    else if (cv2.is_vert)
    {
      // Only cv2 is vertical:
      Comparison_result res1 = compare_y_at_x_2_object()(cv2.ps,cv1.line);
      Comparison_result res2 = compare_y_at_x_2_object()(cv2.pt,cv1.line);

      if (res1 == res2)
      {
        CGAL_assertion(res1 != EQUAL);
        return ((res1 == LARGER) ? SMALLER : LARGER);
      }
      else
        return (EQUAL);
    }

    // Compare using the supporting lines.
    return (compare_y_at_x_2_object()(q, cv1.line, cv2.line));
  }

  /*!
   * Compares the y value of two curves in an epsilon environment to the left
   * of the x-value of their intersection point.
   * \param cv1 The first curve.
   * \param cv2 The second curve.
   * \param q The point.
   * \pre The point q is in the x range of the two curves, and both of them 
   * must be also be defined to its left. Furthermore, cv1(x(q) == cv2(x(q)).
   * \return The relative position of cv1 with respect to cv2 to the left of
   * x(q): LARGER, SMALLER or EQUAL.
   */
  Comparison_result curves_compare_y_at_x_left(const X_monotone_curve_2 & cv1,
                                               const X_monotone_curve_2 & cv2, 
                                               const Point_2 & q) const 
  {
    // The two curve must be defined at q and also to its left.
    CGAL_precondition_code(Compare_xy_2 cmp_xy = compare_xy_2_object(););

    CGAL_precondition((cmp_xy(cv1.ps, q) == SMALLER) ||
                      (cmp_xy(cv1.pt, q) == SMALLER));
    CGAL_precondition((cmp_xy(cv1.ps, q) != SMALLER) ||
                      (cmp_xy(cv1.pt, q) != SMALLER));
    
    CGAL_precondition((cmp_xy(cv2.ps, q) == SMALLER) ||
                      (cmp_xy(cv2.pt, q) == SMALLER));
    CGAL_precondition((cmp_xy(cv2.ps, q) != SMALLER) ||
                      (cmp_xy(cv2.pt, q) != SMALLER));
    
    // Notice q is a placeholder for the x coordinate of the two segments.
    // That is, if we compare them at x(q) the result should be EQUAL.
    CGAL_precondition((is_vertical_2_object()(cv1.line) &&
                       has_on_2_object()(cv1.line, q)) ||
                      (is_vertical_2_object()(cv2.line) &&
                       has_on_2_object()(cv2.line, q)) ||
                      (compare_y_at_x_2_object()(q, cv1.line, cv2.line) ==
                       EQUAL));
    
    // Compare the slopes of the two segments to determine thir relative
    // position immediately to the left of q.
    // Notice we use the supporting lines in order to compare the slopes.
    return (compare_slope_2_object()(cv2.line, cv1.line));
  }

  /*!
   * Compares the y value of two curves in an epsilon environment to the right
   * of the x-value of their intersection point.
   * \param cv1 The first curve.
   * \param cv2 The second curve.
   * \param q The point.
   * \pre The point q is in the x range of the two curves, and both of them 
   * must be also be defined to its right. Furthermore, cv1(x(q) == cv2(x(q)).
   * \return The relative position of cv1 with respect to cv2 to the right of
   * x(q): LARGER, SMALLER or EQUAL.
   */
  Comparison_result
  curves_compare_y_at_x_right(const X_monotone_curve_2 & cv1,
                              const X_monotone_curve_2 & cv2, 
                              const Point_2 & q) const 
  {
    // The two curve must be defined at q and also to its right.
    CGAL_precondition_code(Compare_xy_2 cmp_xy = compare_xy_2_object(););

    CGAL_precondition((cmp_xy(q, cv1.ps) == SMALLER) ||
                      (cmp_xy(q, cv1.pt) == SMALLER));
    CGAL_precondition((cmp_xy(q, cv1.ps) != SMALLER) ||
                      (cmp_xy(q, cv1.pt) != SMALLER));
    
    CGAL_precondition((cmp_xy(q, cv2.ps) == SMALLER) ||
                      (cmp_xy(q, cv2.pt) == SMALLER));
    CGAL_precondition((cmp_xy(q, cv2.ps) != SMALLER) ||
                      (cmp_xy(q, cv2.pt) != SMALLER));

    // Notice q is a placeholder for the x coordinate of the two segments.
    // That is, if we compare them at x(q) the result should be EQUAL.
    CGAL_precondition((is_vertical_2_object()(cv1.line) &&
                       has_on_2_object()(cv1.line, q)) ||
                      (is_vertical_2_object()(cv2.line) &&
                       has_on_2_object()(cv2.line, q)) ||
                      (compare_y_at_x_2_object()(q, cv1.line, cv2.line) ==
                       EQUAL));
    
    // Compare the slopes of the two segments to determine thir relative
    // position immediately to the right of q.
    // Notice we use the supporting lines in order to compare the slopes.
    return (compare_slope_2_object()(cv1.line, cv2.line));
  }
    
  /*! 
   * Return the location of the given point with respect to the input curve.
   * \param cv The curve.
   * \param p The point.
   * \pre p is in the x-range of cv.
   * \return SMALLER if y(p) < cv(x(p)), that is the point is below the curve;
   *         LARGER if y(p) > cv(x(p)), that is the point is above the curve;
   *         or else (if p is on the curve) EQUAL.
   */
  Comparison_result curve_compare_y_at_x(const Point_2 & p,
                                         const X_monotone_curve_2 & cv) const
  {
    CGAL_precondition(point_in_x_range(cv, p));

    if (! cv.is_vert)
    {
      // Compare with the supporting line.
      return compare_y_at_x_2_object()(p, cv.line);
    }
    else
    {
      // Compare with the vertical segment's end-points.
      Compare_y_2       compare_y = compare_y_2_object();
      Comparison_result res1 = compare_y (p, cv.ps);
      Comparison_result res2 = compare_y (p, cv.pt);
      
      if (res1 == res2)
        return (res1);
      else
        return (EQUAL);
    }
  }

  /*! 
   * Check if the two curves are the same (have the same graph).
   * \param cv1 The first curve.
   * \param cv2 The second curve.
   * \return (true) if the two curves are the same.
   */
  bool curve_equal(const X_monotone_curve_2 & cv1,
                   const X_monotone_curve_2 & cv2) const
  {
    Equal_2 equal = equal_2_object();
    return ((equal(cv1.ps, cv2.ps) && equal(cv1.pt, cv2.pt)) ||
            (equal(cv1.ps, cv2.pt) && equal(cv1.pt, cv2.ps)));
  }

  /*!
   * Check if the two points are the same.
   * \param p1 The first point.
   * \param p2 The second point.
   * \return (true) if p1 == p2.
   */
  bool point_equal(const Point_2 & p1, const Point_2 & p2) const
  {
    return (equal_2_object()(p1, p2));
  }
  
  /*!
   * Get the curve source.
   * \param cv The curve.
   * \return The source point.
   */
  const Point_2 & curve_source(const X_monotone_curve_2 & cv) const 
  { 
    return (cv.ps);
  }

  /*!
   * Get the curve target.
   * \param cv The curve.
   * \return The target point.
   */
  const Point_2 & curve_target(const X_monotone_curve_2 & cv) const 
  { 
    return (cv.pt);
  }

  /*!
   * Check whether the curve is x-monotone.
   * \param cv The curves.
   * \return (true) if the curve is x-monotone. In case of segments, the
   * function always returns (true), since all segments are x-monotone. 
   * Vertical segments are also considered as 'weakly' x-monotone.
   */
  bool is_x_monotone(const Curve_2 &) const
  {
    // Return true, since a sgement is always x-monotone.
    return (true);
  }
  
  /*! 
   * Cut the given curve into x-monotone subcurves and insert them to the
   * given output iterator. While segments are x_monotone, still need to pass
   * them out.
   * \param cv The curve.
   * \param o The output iterator
   * \return The past-the-end iterator
   */
  template<class OutputIterator>
  OutputIterator curve_make_x_monotone(const Curve_2 & cv,
                                       OutputIterator o) const
  {
    *o++ = cv;
    return o;
  } 

  /*!
   * Flip a given curve.
   * \param cv The input curve.
   * \return The flipped curve. In case of segments, if the input is [s,t],
   * then the flipped curve is simply [t,s].
   */
  X_monotone_curve_2 curve_opposite(const X_monotone_curve_2 & cv) const
  {
    X_monotone_curve_2 flip_cv(cv);

    flip_cv.ps = cv.pt;
    flip_cv.pt = cv.ps;
    return (flip_cv);
  }

  /*!
   * Split a given curve at a given split point into two sub-curves.
   * \param cv the curve to split
   * \param c1 the output first part of the split curve. Its source is the
   * source of the original curve.
   * \param c2 the output second part of the split curve. Its target is the
   * target of the original curve.
   * \param p the split point.
   * \pre p lies on cv but is not one of its end-points.
   */
  void curve_split(const X_monotone_curve_2 & cv, 
                   X_monotone_curve_2 & c1, X_monotone_curve_2 & c2, 
                   const Point_2 & p) const
  {
    // Check preconditions.
    CGAL_precondition(curve_compare_y_at_x(p, cv) == EQUAL);
    CGAL_precondition_code(Equal_2 is_equal = equal_2_object());
    CGAL_precondition(!is_equal(cv.ps, p));
    CGAL_precondition(!is_equal(cv.pt, p));

    // Do the split.
    c1.line = cv.line;
    c1.ps = cv.ps;
    c1.pt = p;
    c1.is_vert = cv.is_vert;

    c2.line = cv.line;
    c2.ps = p;
    c2.pt = cv.pt;
    c2.is_vert = cv.is_vert;

    return;
  }

  /*!
   * Find the nearest intersection of the two given curves to the right of 
   * a given reference point.
   * Nearest is defined as the lexicographically nearest point, not including 
   * the point reference point itself.
   * If the intersection of the two curves is an X_monotone_curve_2, that is,
   * there is an overlapping subcurve, that contains the reference point in
   * its x-range, the function should return an X_monotone_curve_2 whose 
   * interior is strictly to the right of the reference point (that is, whose
   * left endpoint is the projection of the reference point onto the 
   * overlapping subcurve).
   * \param cv1 The first curve.
   * \param cv2 The second curve.
   * \param p The refernece point.
   * \return An empty object if there is no intersection to the right of p.
   *         An object wrapping a Point_2 in case of a simple intersection.
   *         An object wrapping an X_monotone_curve_2 in case of an overlap.
   */
  Object nearest_intersection_to_right (const X_monotone_curve_2 & cv1,
                                        const X_monotone_curve_2 & cv2,
                                        const Point_2 & p) const
  {
    Point_2  p1, p2;
    bool     is_overlap;

    // Return an empty object if there is no intersection.
    if (! _find_intersection (cv1, cv2, is_overlap, p1, p2))
      return Object();

    // Check if there is a single intersection point.
    if (! is_overlap) 
    {
      // Return the point if it is lexicographically to the right of p.
      if (compare_xy_2_object()(p1, p) == LARGER)
        return (CGAL::make_object (p1));
      
      return Object();
    }
    
    // In case the intersection is an overlapping segment [p1, p2]
    // (notice that p1 < p2):
    if (compare_xy_2_object()(p1, p) == LARGER)
    {
      // The entire segment p1 -> p2 is to the right of p:
      return (CGAL::make_object (X_monotone_curve_2 (cv1.line,
                                                     p1, p2)));
    }
    else if (compare_xy_2_object()(p2, p) == LARGER)
    {
      if (has_on_2_object() (cv1.line, p))
      {
        // p is one the overlapping segment, return it as the first point.
        p1 = p;
      }
      else
      {
        // Perform vertical ray-shooting from p to the overlapping segment
        // and make p1 the resulting point.
        _vertical_ray_shoot (p, cv1,
                             p1);
      }

      // If after the trimming we have p1 == p2, return just a single point.
      if (equal_2_object() (p1, p2))
        return (CGAL::make_object (p1));

      // Return the segment p1 -> p2.
      return (CGAL::make_object (X_monotone_curve_2 (cv1.line,
                                                     p1, p2)));
    }

    // The overlap is entirely to the left of p:
    return Object();
  }

  /*!
   * Find the nearest intersection of the two given curves to the left of 
   * a given reference point.
   * Nearest is defined as the lexicographically nearest point, not including 
   * the point reference point itself.
   * If the intersection of the two curves is an X_monotone_curve_2, that is,
   * there is an overlapping subcurve, that contains the reference point in
   * its x-range, the function should return an X_monotone_curve_2 whose 
   * interior is strictly to the left of the reference point (that is, whose
   * right endpoint is the projection of the reference point onto the 
   * overlapping subcurve).
   * \param cv1 The first curve.
   * \param cv2 The second curve.
   * \param p The refernece point.
   * \return An empty object if there is no intersection to the left of p.
   *         An object wrapping a Point_2 in case of a simple intersection.
   *         An object wrapping an X_monotone_curve_2 in case of an overlap.
   */
  Object nearest_intersection_to_left (const X_monotone_curve_2 & cv1,
                                       const X_monotone_curve_2 & cv2,
                                       const Point_2 & p) const
  {
    Point_2  p1, p2;
    bool     is_overlap;

    // Return an empty object if there is no intersection.
    if (! _find_intersection (cv1, cv2, is_overlap, p1, p2))
      return Object();

    // Check if there is a single intersection point.
    if (! is_overlap) 
    {
      // Return the point if it is lexicographically to the left of p.
      if (compare_xy_2_object()(p1, p) == SMALLER)
        return (CGAL::make_object (p1));

      return Object();
    }

    // In case the intersection is an overlapping segment [p1, p2]
    // (notice that p1 < p2):
    if (compare_xy_2_object()(p2, p) == SMALLER)
    {
      // The entire segment p1 -> p2 is to the left of p:
      return (CGAL::make_object (X_monotone_curve_2 (cv1.line,
                                                     p1, p2)));
    }
    else if (compare_xy_2_object()(p1, p) == SMALLER)
    {
      if (has_on_2_object() (cv1.line, p))
      {
        // p is one the overlapping segment, return it as the first point.
        p2 = p;
      }
      else
      {
        // Perform vertical ray-shooting from p to the overlapping segment
        // and make p2 the resulting point.
        _vertical_ray_shoot (p, cv1,
                             p2);
      }

      // If after the trimming we have p1 == p2, return just a single point.
      if (equal_2_object() (p1, p2))
        return (CGAL::make_object (p1));

      // Return the segment p1 -> p2.
      return (CGAL::make_object (X_monotone_curve_2 (cv1.line,
                                                     p1, p2)));
    }
     
    // The overlap is entirely to the right of p:
    return Object();
  }

  /*!
   * Check whether the two given curves overlap.
   * \patam cv1 The first curve.
   * \patam cv2 The second curve.
   * \return (true) if the two curves overlap in a one-dimensional subcurve
   * (i.e., not in a finite number of points). Otherwise, if they have a finite
   * number of intersection points, or none at all, return (false).
   */
  bool curves_overlap(const X_monotone_curve_2 & cv1,
                      const X_monotone_curve_2 & cv2) const
  {
    // Comparing lines
    //   if (!equal_2_object()(cv1.line, cv2.line)) return false;
    // doesn't work, as coincident lines with opposite direction are
    // considered different!
    Has_on_2 has_on = has_on_2_object();
    
    if (!has_on(cv1.line, cv2.ps) || !has_on(cv1.line, cv2.pt))
      return false;

    if (cv1.is_vert)
    {
      if (cv2.is_vert)
      {
        Compare_y_2 compare_y = compare_y_2_object();
        Comparison_result res_ss = compare_y (cv1.ps, cv2.ps);
        Comparison_result res_st = compare_y (cv1.ps, cv2.pt);

        if (res_ss == SMALLER)
        {
          if (res_st == LARGER)
            return true;
          
          if (compare_y (cv1.pt, cv2.ps) == LARGER)
            return true;
          
          return (compare_y (cv1.pt, cv2.pt) == LARGER);
        }

        if (res_ss == LARGER)
        {
          if (res_st == SMALLER)
            return true;
          
          if (compare_y (cv1.pt, cv2.ps) == SMALLER)
            return true;
          
          return (compare_y (cv1.pt, cv2.pt) == SMALLER);
        }

        // res_ss == EQUAL
        if (res_st == SMALLER)
          return (compare_y (cv1.pt, cv2.ps) == LARGER);
        
        return (compare_y (cv1.pt, cv2.ps) == SMALLER);
      }
      return false;
    }

    if (cv2.is_vert)
      return false;

    Compare_x_2 compare_x = compare_x_2_object();
    Comparison_result res_ss = compare_x (cv1.ps, cv2.ps);
    Comparison_result res_st = compare_x (cv1.ps, cv2.pt);
    
    if (res_ss == SMALLER)
    {
      if (res_st == LARGER)
        return true;
      
      if (compare_x (cv1.pt, cv2.ps) == LARGER)
        return true;
      
      return (compare_x (cv1.pt, cv2.pt) == LARGER);
    }

    if (res_ss == LARGER)
    {
      if (res_st == SMALLER)
        return true;
      
      if (compare_x (cv1.pt, cv2.ps) == SMALLER)
        return true;
      
      return (compare_x (cv1.pt, cv2.pt) == SMALLER);
    }

    // res_ss == EQUAL
    if (res_st == SMALLER)
      return (compare_x (cv1.pt, cv2.ps) == LARGER);
    
    return (compare_x (cv1.pt, cv2.ps) == SMALLER);
  }

private:

  /*!
   * Find the intersection between teo segements.
   * \param cv1 The first segment.
   * \param cv2 The second segment.
   * \param is_ovelap Are the two segment overlapping.
   * \param p1 The intersection point (if there is no overlap),
   * otherwise the leftmost end-point of the intersection segment.
   * \param p2 If there is an overlap, the rightmost end-point of the 
   * intersection segment.
   * \return (true) if an intersection has been found.
   */
  bool _find_intersection (const X_monotone_curve_2 & cv1,
                           const X_monotone_curve_2 & cv2,
                           bool & is_overlap,
                           Point_2 & p1, Point_2 & p2) const
  {
    // Computing the orientation ahead and checking whether the end points of
    // one curve are in opposite orientations with respect to the other seems
    // slow down the process!

    // Intersect the two supporting lines.
    Object    res = intersect_2_object()(cv1.line, cv2.line);

    // Parallel lines do not intersect
    if (res.is_empty())
      return (false);
    
    Point_2   ip;
    if (assign(ip, res))
    {
      is_overlap = false;

      // Simple case of intersection at a single point.
      if (_is_on_segment(cv1, ip) && _is_on_segment(cv2, ip))
      {
        p1 = ip;
        return (true);
      }
      return (false);
    }
    else // The two supporting lines overlap. 
    {
      // Assign the end-points such that p1 < p2.
      Compare_xy_2       comp_xy = compare_xy_2_object();

      if (comp_xy(cv1.ps, cv1.pt) == SMALLER)
      {
        p1 = cv1.ps;
        p2 = cv1.pt;
      }
      else
      {
        p1 = cv1.pt;
        p2 = cv1.ps;
      }

      // Clip the first segment with respect to cv2.
      Comparison_result res = comp_xy(cv2.ps, cv2.pt);
      const Point_2 &   left2 = (res == SMALLER) ? cv2.ps : cv2.pt;
      const Point_2 &   right2 = (res == LARGER) ? cv2.ps : cv2.pt;

      if (comp_xy(p2, left2) == SMALLER)
        return (false);
      else if (comp_xy(p1, left2) == SMALLER)
        p1 = left2;

      if (comp_xy(p1, right2) == LARGER)
        return (false);
      else if (comp_xy(p2, right2) == LARGER)
        p2 = right2;

      // Check if the intersection segment has not become a point.
      is_overlap = (comp_xy(p1,p2) != EQUAL);
      CGAL_assertion(comp_xy(p1,p2) != LARGER);
      return (true);
    }
  }

  /*!
   * Check whether the given point lies on the given segment.
   * \param cv The curve.
   * \param q The point.
   * \pre The function assumes q lies on the supporting line of the segment.
   * \return (true) if q lies of cv.
   */
  bool _is_on_segment (const X_monotone_curve_2 & cv, const Point_2 & q) const
  {
    if (! cv.is_vert)
    {
      Compare_x_2       comp_x = compare_x_2_object();
      Comparison_result res1 = comp_x (q, cv.ps);
      Comparison_result res2 = comp_x (q, cv.pt);

      // We check if x(q) equals the x value of one of the end-points.
      // If not, we check whether one end-point is to q's left and the other is
      // to its right.
      return ((res1 == EQUAL) || (res2 == EQUAL) ||
              (res1 != res2));
    }
    else
    {
      Compare_y_2       comp_y = compare_y_2_object();
      Comparison_result res1 = comp_y (q, cv.ps);
      Comparison_result res2 = comp_y (q, cv.pt);

      // We check if x(q) equals the y value of one of the end-points.
      // If not, we check whether one end-point is above q and the other is
      // below it.
      return ((res1 == EQUAL) || (res2 == EQUAL) ||
              (res1 != res2));
    }
  }

  /*!
   * Perform vertical ray-shooting from a given point towards a given curve.
   * \param q The source point of the ray.
   * \param cv The target curve.
   * \param p The resulting point.
   * \return Whether we have successfully computed a point p.
   */
  bool _vertical_ray_shoot (const Point_2& q, const X_monotone_curve_2& cv,
                            Point_2& p) const
  {
    // Construct a vertical line passing through q.
    typename Kernel::Direction_2  dir (0, 1);
    typename Kernel::Line_2       vl = construct_line_2_object() (q, dir);

    // Compute the intersetion between the vertical line and the line
    // supporting the curve cv.
    Object    res = intersect_2_object()(cv.line, vl);
    bool      ray_shoot_successful = assign(p, res);

    CGAL_assertion (ray_shoot_successful);
    return (ray_shoot_successful);
  }
};

/*!
 * A representation of a segment, as used by the Arr_segment_cached_traits_2
 * traits class.
 */
template <class Kernel_>
class Segment_cached_2 :
    public Arr_segment_cached_traits_2<Kernel_>::My_segment_cached_2
{
  typedef typename Arr_segment_cached_traits_2<Kernel_>::My_segment_cached_2
                                                  Base;
  typedef typename Kernel_::Segment_2             Segment_2;
  typedef typename Kernel_::Point_2               Point_2;
  typedef typename Kernel_::Line_2                Line_2;

#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_3
public:
  using Base::ps;
  using Base::pt;
#endif

  public:

  /*!
   * Default constructor.
   */
  Segment_cached_2 () :
    Base()
  {}
    
  /*!
   * Constructor from a segment.
   * \param seg The segment.
   */
  Segment_cached_2 (const Segment_2 & seg) :
    Base(seg)
  {}

  /*!
   * Construct a segment from two end-points.
   * \param source The source point.
   * \param target The target point.
   */
  Segment_cached_2 (const Point_2 & source, const Point_2 & target) :
    Base(source,target)
  {}

  /*!
   * Construct a segment from a line and two end-points.
   * \param line The supporting line.
   * \param source The source point.
   * \param target The target point.
   * \pre Both source and target must be on the supporting line.
   */
  Segment_cached_2 (const Line_2 & line,
                    const Point_2 & source, const Point_2 & target) :
    Base(line,source,target)
  {}

  /*!
   * Cast to a segment.
   */
  operator Segment_2 () const
  {
    return (Segment_2(ps, pt));
  }

  /*!
   * Create a bounding box for the segment.
   */
  Bbox_2 bbox() const
  {
    Segment_2 seg(ps, pt);
    return (seg.bbox());
  }

  /*!
   * Get the segment source.
   */
  const Point_2 & source() const 
  { 
    return ps; 
  }

  /*!
   * Get the segment target.
   */
  const Point_2 & target() const
  { 
    return pt;
  }
};

/*!
 * Exporter for a cached segment.
 */
template <class Kernel_, class Stream_>
Stream_ & operator<<(Stream_ & os, const Segment_cached_2<Kernel_> & seg)
{
  os << static_cast<typename Kernel_::Segment_2>(seg);
  return (os);
}

/*!
 * Importer for a cached segment.
 */
template <class Kernel_, class Stream_>
Stream_ & operator>>(Stream_ & is, Segment_cached_2<Kernel_> & seg)
{
  typename Kernel_::Segment_2 kernel_seg;
  is >> kernel_seg;
  seg = kernel_seg;
  return is;
}

CGAL_END_NAMESPACE

#endif
