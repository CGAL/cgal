// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $$
// release_date  : $$
//
// file          : include/CGAL/Arr_segment_cached_traits_2.h
// package       : Planar_map (5.87)
// maintainer    : Efi Fogel         <efif@math.tau.ac.il>
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Ron Wein          <wein@post.tau.ac.il>
//                 Efi Fogel         <efif@post.tau.ac.il>
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ======================================================================
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
    My_segment_cached_2 (const Segment_2& seg)
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
    My_segment_cached_2 (const Point_2& source, const Point_2& target)
    {
      Kernel_   kernel;
      
      line = kernel.construct_line_2_object()(source, target);
      is_vert = kernel.is_vertical_2_object()(line);
      
      ps = source;
      pt = target;
    }

    friend class Arr_segment_cached_traits_2;
  };

  // Traits objects
  typedef typename Kernel::Point_2        Point_2;
  typedef Segment_cached_2<Kernel>        X_curve_2;
  typedef Segment_cached_2<Kernel>        Curve_2;

  // Backward compatability    
  typedef Point_2                         Point;
  typedef X_curve_2                       X_curve;
  typedef X_curve_2                       Curve;

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
  bool curve_is_vertical(const X_curve_2 & cv) const 
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
  bool curve_is_in_x_range(const X_curve_2 & cv, const Point_2 & q) const
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
  Comparison_result curve_compare_at_x(const X_curve_2 & cv1, 
				       const X_curve_2 & cv2, 
				       const Point_2 & q) const
  {
    CGAL_precondition(curve_is_in_x_range(cv1, q));
    CGAL_precondition(curve_is_in_x_range(cv2, q));

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
	const Point_2&    lower1 = (res1 == SMALLER) ? cv1.ps : cv1.pt;
	const Point_2&    upper1 = (res1 == SMALLER) ? cv1.pt : cv1.ps;
	Comparison_result res2 = compare_y (cv2.ps, cv2.pt);
	const Point_2&    lower2 = (res2 == SMALLER) ? cv2.ps : cv2.pt;
	const Point_2&    upper2 = (res2 == SMALLER) ? cv2.pt : cv2.ps;

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
  Comparison_result curve_compare_at_x_left(const X_curve_2 & cv1,
                                            const X_curve_2 & cv2, 
                                            const Point_2 & q) const 
  {
    // The two curves must not be vertical.
    CGAL_precondition(! cv1.is_vert);
    CGAL_precondition(! cv2.is_vert);

    // The two curve must be defined at q and also to its left.
    CGAL_precondition_code(
	Less_x_2 less_x = less_x_2_object();
	);

    CGAL_precondition (less_x(cv1.ps, q) || less_x(cv1.pt, q));
    CGAL_precondition (!(less_x(cv1.ps, q) && less_x(cv1.pt, q)));
    
    CGAL_precondition (less_x(cv2.ps, q) || less_x(cv2.pt, q));
    CGAL_precondition (!(less_x(cv2.ps, q) && less_x(cv2.pt, q)));
    
    // Notice q is a placeholder for the x coordinate of the two segments.
    // That is, if we compare them at x(q) the result should be EQUAL.
    CGAL_precondition(compare_y_at_x_2_object()(q,cv1.line,cv2.line) == EQUAL);
    
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
  Comparison_result curve_compare_at_x_right(const X_curve_2 & cv1,
                                             const X_curve_2 & cv2, 
                                             const Point_2 & q) const 
  {
    // The two curves must not be vertical.
    CGAL_precondition(! cv1.is_vert);
    CGAL_precondition(! cv2.is_vert);

    // The two curve must be defined at q and also to its right.
    CGAL_precondition_code(Less_x_2 less_x = less_x_2_object(););

    CGAL_precondition (less_x(q, cv1.ps) || less_x(q, cv1.pt));
    CGAL_precondition (!(less_x(q, cv1.ps) && less_x(q, cv1.pt)));
    
    CGAL_precondition (less_x(q, cv2.ps) || less_x(q, cv2.pt));
    CGAL_precondition (!(less_x(q, cv2.ps) && less_x(q, cv2.pt)));

    // Notice q is a placeholder for the x coordinate of the two segments.
    // That is, if we compare them at x(q) the result should be EQUAL.
    CGAL_precondition(compare_y_at_x_2_object()(q,cv1.line,cv2.line) == EQUAL);
    
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
   * \return SMALLER if cv(x(p)) < y(p);
   *         LARGER if cv(x(p)) > y(p);
   *         or else (if p is on the curve) EQUAL.
   */
  Comparison_result curve_get_point_status (const X_curve_2 & cv, 
					    const Point_2 & p) const
  {
    CGAL_precondition(curve_is_in_x_range(cv, p));

    if (! cv.is_vert)
    {
      // Compare with the supporting line.
      Comparison_result res = compare_y_at_x_2_object()(p, cv.line);

      if (res == LARGER)
	return (SMALLER);
      else if (res == SMALLER)
	return (LARGER);
      return (EQUAL);
    }
    else
    {
      // Compare with the vertical segment's end-points.
      Compare_y_2       compare_y = compare_y_2_object();
      Comparison_result res1 = compare_y (cv.ps, p);
      Comparison_result res2 = compare_y (cv.pt, p);
      
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
  bool curve_is_same(const X_curve_2 & cv1,const X_curve_2 & cv2) const
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
  bool point_is_same(const Point_2 & p1, const Point_2 & p2) const
  {
    return (equal_2_object()(p1, p2));
  }
  
  /*!
   * Get the curve source.
   * \param cv The curve.
   * \return The source point.
   */
  const Point_2& curve_source(const X_curve_2 & cv) const 
  { 
    return (cv.ps);
  }

  /*!
   * Get the curve target.
   * \param cv The curve.
   * \return The target point.
   */
  const Point_2& curve_target(const X_curve_2 & cv) const 
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
  
  /*! Cut the given curve into x-monotone subcurves and stores them in the
   * given list. While segments are x_monotone, still need to cast their type.
   * \param cv The curve.
   * \param x_curves A list of the output x-monotone sub-curves.
   */
  void make_x_monotone(const Curve_2& cv, std::list<Curve_2>& l) const
  {
    l.clear();
    l.push_back(X_curve_2(cv));
  } 

  /*!
   * Flip a given curve.
   * \param cv The input curve.
   * \return The flipped curve. In case of segments, if the input is [s,t],
   * then the flipped curve is simply [t,s].
   */
  X_curve_2 curve_flip(const X_curve_2 & cv) const
  {
    X_curve_2 flip_cv(cv);

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
  void curve_split(const X_curve_2& cv, 
		   X_curve_2& c1, X_curve_2& c2, 
                   const Point_2& p) const
  {
    // Check preconditions.
    CGAL_precondition(curve_get_point_status(cv, p) == EQUAL);
    CGAL_precondition_code(Compare_xy_2 compare_xy = compare_xy_2_object());
    CGAL_precondition(compare_xy(cv.ps, p) != EQUAL);
    CGAL_precondition(compare_xy(cv.pt, p) != EQUAL);
    
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
   * Find the nearest intersection point (or points) of two given curves to
   * the right lexicographically of a given point not includin the point
   * itself, (with one exception explained below).
   * If the intersection of the two curves is an X_curve_2, that is, they
   * overlap at infinitely many points, then if the right endpoint and the
   * left endpoint of the overlapping subcurve are strickly to the right of
   * the given point, they are returned through the two other point
   * references respectively. If the given point is between the
   * overlapping-subcurve endpoints, or the point is its left endpoint,
   * the point and the right endpoint of the subcurve are returned through
   * the point references respectively. If the intersection of the two curves
   * is a point to the right of the given point, it is returned through the
   * point references.
   * \param cv1 The first curve.
   * \param cv2 The second curve.
   * \param p The refernece point.
   * \param p1 The first output point.
   * \param p2 The second output point.
   * \return (true) if c1 and c2 do intersect to the right of p, or (false)
   * if no such intersection exists.
   */
  bool nearest_intersection_to_right(const X_curve_2 & cv1,
                                     const X_curve_2 & cv2,
                                     const Point_2 & p,
                                     Point_2 & p1, Point_2 & p2) const
  {
    bool     is_overlap;

    if (! _find_intersection (cv1, cv2, is_overlap, p1, p2))
      return (false);

    if (! is_overlap) 
    {
      if (compare_xy_2_object()(p1, p) == LARGER)
      {
	p2 = p1;
	return (true);
      }
      return (false);
    }
    else
    {
      // Notice that p1 < p2.
      if (compare_xy_2_object()(p1, p) == LARGER)
      {
	return (true);
      }
      else if (compare_xy_2_object()(p2, p) == LARGER)
      {
	p1 = p;
	return (true);
      }
      return (false);
    }
  }

  /*!
   * Find the nearest intersection point of two given curves to the left of 
   * a given point. Nearest is defined as the lexicographically nearest not 
   * including the point itself (with one exception explained below).
   * If the intersection of the two curves is an X_curve_2, that is,
   * there is an overlapping subcurve, then if the the source and target of the
   * subcurve are strickly to the left, they are returned through two
   * other point references p1 and p2. If p is between the source and target
   * of the overlapping subcurve, or p is its right endpoint, p and the source
   * of the left endpoint of the subcurve are returned through p1 and p2 
   * respectively.
   * If the intersection of the two curves is a point to the left of p, it is
   * returned through the p1 and p2.
   * \param cv1 The first curve.
   * \param cv2 The second curve.
   * \param p The refernece point.
   * \param p1 The first output point.
   * \param p2 The second output point.
   * \return (true) if c1 and c2 do intersect to the left of p, or (false)
   * if no such intersection exists.
   */
  bool nearest_intersection_to_left(const X_curve_2 & cv1,
                                    const X_curve_2 & cv2,
                                    const Point_2 & p,
                                    Point_2 & p1, Point_2 & p2) const
  {
    bool     is_overlap;

    if (! _find_intersection (cv1, cv2, is_overlap, p1, p2))
      return (false);

    if (! is_overlap) 
    {
      if (compare_xy_2_object()(p1, p) == SMALLER)
      {
	p2 = p1;
	return (true);
      }
      return (false);
    }
    else
    {
      // Notice that ip1 < ip2.
      if (compare_xy_2_object()(p2, p) == SMALLER)
      {
	return (true);
      }
      else if (compare_xy_2_object()(p1, p) == SMALLER)
      {
	p2 = p;
	return (true);
      }
      return (false);
    }
  }

  /*!
   * Check whether the two given curves overlap.
   * \patam cv1 The first curve.
   * \patam cv2 The second curve.
   * \return (true) if the two curves overlap in a one-dimensional subcurve
   * (i.e., not in a finite number of points). Otherwise, if they have a finite
   * number of intersection points, or none at all, return (false).
   */
  bool curves_overlap(const X_curve_2 & cv1, const X_curve_2 & cv2) const
  {
    // Comparing lines
    //   if (!equal_2_object()(cv1.line, cv2.line)) return false;
    // doesn't work, as coincident lines with opposite direction are
    // considered different!
    Has_on_2 has_on = has_on_2_object();
    if (!has_on(cv1.line, cv2.ps) || !has_on(cv1.line, cv2.pt)) return false;

    if (cv1.is_vert) {
      if (cv2.is_vert) {
        Compare_y_2 compare_y = compare_y_2_object();
        Comparison_result res_ss = compare_y (cv1.ps, cv2.ps);
        Comparison_result res_st = compare_y (cv1.ps, cv2.pt);
        if (res_ss == SMALLER) {
          if (res_st == LARGER) return true;
          if (compare_y (cv1.pt, cv2.ps) == LARGER) return true;
          return (compare_y (cv1.pt, cv2.pt) == LARGER);
        }

        if (res_ss == LARGER) {
          if (res_st == SMALLER) return true;
          if (compare_y (cv1.pt, cv2.ps) == SMALLER) return true;
          return (compare_y (cv1.pt, cv2.pt) == SMALLER);
        }

        // res_ss == EQUAL
        if (res_st == SMALLER)
          return (compare_y (cv1.pt, cv2.ps) == LARGER);
        return (compare_y (cv1.pt, cv2.ps) == SMALLER);
      }
      return false;
    }
    if (cv2.is_vert) return false;

    Compare_x_2 compare_x = compare_x_2_object();
    Comparison_result res_ss = compare_x (cv1.ps, cv2.ps);
    Comparison_result res_st = compare_x (cv1.ps, cv2.pt);
    if (res_ss == SMALLER) {
      if (res_st == LARGER) return true;
      if (compare_x (cv1.pt, cv2.ps) == LARGER) return true;
      return (compare_x (cv1.pt, cv2.pt) == LARGER);
    }

    if (res_ss == LARGER) {
      if (res_st == SMALLER) return true;
      if (compare_x (cv1.pt, cv2.ps) == SMALLER) return true;
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
   * \return (true) if an intersectio has been found.
   */
  bool _find_intersection (const X_curve_2& cv1,
			   const X_curve_2& cv2,
                           bool& is_overlap,
			   Point_2& p1, Point_2& p2) const
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
      const Point_2&    left2 = (res == SMALLER) ? cv2.ps : cv2.pt;
      const Point_2&    right2 = (res == LARGER) ? cv2.ps : cv2.pt;

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
  bool _is_on_segment (const X_curve_2 & cv, const Point_2 & q) const
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
  Segment_cached_2 (const Segment_2& seg) :
    Base(seg)
  {}

  /*!
   * Construct a segment from two end-points.
   * \param source The source point.
   * \param target The target point.
   */
  Segment_cached_2 (const Point_2& source, const Point_2& target) :
    Base(source,target)
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
  Bbox_2 bbox()
  {
    Segment_2 seg(ps, pt);
    return (seg.bbox());
  }

  /*!
   * Get the segment source.
   */
  const Point_2& source() const 
  { 
    return ps; 
  }

  /*!
   * Get the segment target.
   */
  const Point_2& target() const
  { 
    return pt;
  }
};

/*!
 * Output operator for a cached segment.
 */
template <class Kernel_>
::std::ostream& operator<< (::std::ostream& os,
                            const Segment_cached_2<Kernel_>& seg)
{
  os << static_cast<typename Kernel_::Segment_2>(seg);
  return (os);
}

CGAL_END_NAMESPACE

#endif
