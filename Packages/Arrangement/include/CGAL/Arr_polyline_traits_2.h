// Copyright (c) 2003  Tel-Aviv University (Israel).
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
// Author(s)     : Ron Wein <wein@post.tau.ac.il>
#ifndef CGAL_ARR_POLYLINE_TRAITS_H
#define CGAL_ARR_POLYLINE_TRAITS_H

#include <CGAL/basic.h>
#include <CGAL/tags.h>
#include <CGAL/intersections.h>

#include <list>
#include <fstream>
#include <vector>

CGAL_BEGIN_NAMESPACE

template <class Segment_traits_> class Polyline_2;

template <class Segment_traits_>
class Arr_polyline_traits_2
{
public:

  typedef Tag_true                                  Has_left_category;
  typedef Tag_false                                 Has_reflect_category;
    
  typedef Segment_traits_                           Segment_traits_2;
  typedef typename Segment_traits_2::Kernel         Kernel;
  typedef Arr_polyline_traits_2<Segment_traits_2>   Self;
  typedef typename Segment_traits_2::Point_2        Point_2;
  typedef typename Segment_traits_2::Curve_2        Segment_2;

  typedef Polyline_2<Segment_traits_2>              Curve_2;
  typedef Polyline_2<Segment_traits_2>              X_monotone_curve_2;

  // Obsolete, for backward compatibility
  typedef Point_2                       Point;
  typedef X_monotone_curve_2            X_curve;
  typedef Curve_2                       Curve;

protected:

  Segment_traits_2   seg_traits;   // A traits class for segments.

public:

  /*!
   * Default constructor.
   */
  Arr_polyline_traits_2 () :
    seg_traits()
  {}

  /*!
   * Compare the x-coordinates of two given points.
   * \param p1 The first point.
   * \param p2 The second point.
   * \return LARGER if x(p1) > x(p2); SMALLER if x(p1) < x(p2); or else EQUAL.
   */
  Comparison_result compare_x (const Point_2& p1, const Point_2& p2) const
  {
    return (seg_traits.compare_x (p1, p2));
  }

  /*! 
   * Compares lexigoraphically the two points: by x, then by y.
   * \param p1 Te first point.
   * \param p2 The second point.
   * \return LARGER if x(p1) > x(p2), or if x(p1) = x(p2) and y(p1) > y(p2); 
   *         SMALLER if x(p1) < x(p2), or if x(p1) = x(p2) and y(p1) < y(p2);
   *         or else EQUAL.
   */
  Comparison_result compare_xy(const Point_2& p1, const Point_2& p2) const
  {
    return (seg_traits.compare_xy (p1, p2));
  }

  /*!
   * Check whether the given curve is a vertical segment.
   * \param cv The curve.
   * \return (true) if the curve is vertical.
   */
  bool curve_is_vertical(const X_monotone_curve_2& cv) const
  {
    // An x-monotone polyline can represent a vertical segment only if it
    // is comprised of a single segment.
    return (cv._size() == 1 && 
	    seg_traits.curve_is_vertical(cv[0]));
  } 

  /*!
   * Check whether the given point is in the x-range of the given curve.
   * In out case, the curve is a segment [s, t], check whether x(s)<=x(q)<=x(t)
   * or whether x(t)<=x(q)<=x(s).
   * \param cv The curve.
   * \param q The point.
   * \return (true) if q is in the x-range of cv.
   */
  bool point_in_x_range(const X_monotone_curve_2& cv, const Point_2& q) const
  {
    return (_locate_point (cv, q) >= 0);
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
    // Get the indices of the segments in cv1 and cv2 containing q.
    int    i1 = _locate_point (cv1, q);
    int    i2 = _locate_point (cv2, q);

    CGAL_precondition(i1 >= 0);
    CGAL_precondition(i2 >= 0);

    return (seg_traits.curves_compare_y_at_x (cv1[i1], cv2[i2], q));
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
  Comparison_result curves_compare_y_at_x_left(const X_monotone_curve_2& cv1,
                                               const X_monotone_curve_2& cv2, 
                                               const Point_2& q) const 
  {
    // The two curves must not be vertical.
    CGAL_precondition(! curve_is_vertical(cv1));
    CGAL_precondition(! curve_is_vertical(cv2));

    // Get the indices of the segments in cv1 and cv2 containing q and
    // defined to its left.
    int    i1 = _locate_point_side (cv1, q, false);
    int    i2 = _locate_point_side (cv2, q, false);

    CGAL_precondition(i1 >= 0);
    CGAL_precondition(i2 >= 0);

    // Compare cv1[i1] and cv2[i2] at q's left.
    return (seg_traits.curves_compare_y_at_x_left (cv1[i1], cv2[i2], q));
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
  Comparison_result curves_compare_y_at_x_right(const X_monotone_curve_2& cv1,
						const X_monotone_curve_2& cv2, 
						const Point_2& q) const 
  {
    // The two curves must not be vertical.
    CGAL_precondition(! curve_is_vertical(cv1));
    CGAL_precondition(! curve_is_vertical(cv2));

    // Get the indices of the segments in cv1 and cv2 containing q and
    // defined to its right.
    int    i1 = _locate_point_side (cv1, q, true);
    int    i2 = _locate_point_side (cv2, q, true);

    CGAL_precondition(i1 >= 0);
    CGAL_precondition(i2 >= 0);

    // Compare cv1[i1] and cv2[i2] at q's right.
    return (seg_traits.curves_compare_y_at_x_right (cv1[i1], cv2[i2], q));
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
  Comparison_result curve_compare_y_at_x (const Point_2& p,
                                          const X_monotone_curve_2 &cv) const
  {
    // Get the indices of the segments in cv containing p.
    int    i = _locate_point (cv, p);

    CGAL_precondition(i >= 0);

    // Compare the segment cv[i] and p.
    return (seg_traits.curve_compare_y_at_x (p, cv[i]));
  }
    
  /*!
   * Check if the two points are the same.
   * \param p1 The first point.
   * \param p2 The second point.
   * \return (true) if p1 == p2.
   */
  bool point_equal (const Point_2& p1, const Point_2& p2) const
  { 
    return (seg_traits.point_equal(p1, p2));
  }
  
  /*! 
   * Check if the two curves are the same (have the same graph).
   * \param cv1 The first curve.
   * \param cv2 The second curve.
   * \return (true) if the two curves are the same.
   */
  bool curve_equal(const X_monotone_curve_2& cv1,
                   const X_monotone_curve_2& cv2) const
  {
    // The two curves must contain the same number of segments.
    int    n1 = cv1._size();
    int    n2 = cv2._size();

    if (n1 != n2) 
      return (false);

    // Check the pairwise equality of the contained segments. 
    bool   equal = true;
    int    i;

    for (i = 0; i < n1; i++)
    {
      if (! seg_traits.curve_equal (cv1[i], cv2[i]))
      {
	equal = false;
	break;
      }
    }

    if (equal)
      return (true);

    // Check the reverse order.
    for (i = 0; i < n1; i++)
    {
      if (! seg_traits.curve_equal (cv1[i], cv2[n1 - i - 1]))
	return (false);
    }

    return (true);
  }

  /*!
   * Get the curve source (the source of the first segment in the polyline).
   * \param cv The curve.
   * \return The source point.
   */
  Point_2 curve_source (const X_monotone_curve_2& cv) const
  {
    CGAL_assertion(cv._size() > 0);
    return (seg_traits.curve_source(cv[0]));
  }

  /*!
   * Get the curve target (the target of the last segment in the polyline).
   * \param cv The curve.
   * \return The target point.
   */
  Point_2 curve_target (const X_monotone_curve_2 & cv) const
  {
    CGAL_assertion(cv._size() > 0);
    return (seg_traits.curve_target(cv[cv._size() - 1]));
  }

  /*!
   * Flip a given curve.
   * \param cv The input curve.
   * \return The flipped curve.
   */
  X_monotone_curve_2 curve_opposite (const X_monotone_curve_2& cv) const
  {
    // Copy the segments in revered order while flipping each one individually.
    int                n = cv._size();
    int                i;
    X_monotone_curve_2 flip_cv;

    for (i = 0; i < n; i++)
    {
      flip_cv._push_back (seg_traits.curve_opposite(cv[n - i - 1]));
    }

    return (flip_cv); 
  }

  /*!
   * Check whether the curve is x-monotone.
   * \param curve The curve.
   * \return (true) if the curve is x-monotone.
   */
  bool is_x_monotone(const Curve_2& curve) const
  {
    // Go over all curve segments.
    int               n = curve._size();
    Comparison_result curr_res;
    Comparison_result prev_res = EQUAL;
    int               i;

    for (i = 0; i < n; i++)
    {
      curr_res = seg_traits.compare_x (seg_traits.curve_source(curve[i]),
				       seg_traits.curve_target(curve[i]));

      // In case of a vertical segment, it must be the only segment in the 
      // polyline to be considered x-monotone.
      if (curr_res == EQUAL)
	return (curve._size() == 1);

      if (curr_res != prev_res && prev_res != EQUAL)
	return (false);
      
      prev_res = curr_res;
    }

    return (true);
  }
  
  /*! 
   * Cut the given curve into x-monotone subcurves and insert them to the
   * given output iterator.
   * \param cv The curve.
   * \param o The output iterator
   * \return The past-the-end iterator
   */
  template<class OutputIterator>
  OutputIterator curve_make_x_monotone (const Curve_2& curve,
					OutputIterator o) const
  {
    // Go over all curve segments.
    int                n = curve._size();
    Comparison_result  curr_res;
    Comparison_result  prev_res = EQUAL;
    X_monotone_curve_2 x_cv;
    int                i;

    for (i = 0; i < n; i++)
    {
      curr_res = seg_traits.compare_x (seg_traits.curve_source(curve[i]),
				       seg_traits.curve_target(curve[i]));

      // In case of a vertical segment, it must be the only segment in the 
      // polyline to be considered x-monotone.
      if (curr_res == EQUAL)
      {
	// Current segment is vertical:
	if (x_cv._size() != 0)
	{
	  // Cut the previous chain.
	  *o++ = x_cv;
	  x_cv._clear();
	}

	// Insert the vertical segment as a singleton polyline.
	x_cv._push_back (curve[i]);
	*o++ = x_cv;
	x_cv._clear();
      }
      else if (curr_res != prev_res && prev_res != EQUAL)
      {
	// We should cut the chain at this point and start a new one.
	if (x_cv._size() != 0)
	{
	  // Cut the previous chain.
	  *o++ = x_cv;
	  x_cv._clear();
	}
	
	x_cv._push_back (curve[i]);
      }
      else
      {
	// Just append [ps, pt] to the current chain.
	x_cv._push_back (curve[i]);
      }
      
      prev_res = curr_res;
    }

    // Append the last x-monotone chain.
    if (x_cv._size() != 0)
    {
      // Cut the previous chain.
      *o++ = x_cv;
      x_cv._clear();
    }

    return (o);
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
  void curve_split (const X_monotone_curve_2& cv, 
		    X_monotone_curve_2& c1, X_monotone_curve_2& c2, 
		    const Point_2& p) const
  {
    // Locate the segment on the polyline cv that contains p.
    int    i = _locate_point (cv, p);

    // Check preconditions.
    CGAL_precondition(i >= 0);
    CGAL_precondition(! seg_traits.point_equal(curve_source(cv), p));
    CGAL_precondition(! seg_traits.point_equal(curve_target(cv), p));

    // Clear the output curves.
    c1._clear(); 
    c2._clear();

    // Push all segments labeled (0, 1, ... , i-1) into c1.
    int    j;

    for (j = 0; j < i; j++)
      c1._push_back (cv[j]);

    // Check whether the split point is cv[i]'s source of target.
    if (seg_traits.point_equal (seg_traits.curve_source(cv[i]), p))
    {
      // cv[i] should be the first segment in c2.
      c2._push_back (cv[i]);
    }
    else if (seg_traits.point_equal (seg_traits.curve_target(cv[i]), p))
    {
      // cv[i] should be the last segment in c1.
      c1._push_back (cv[i]);
    }
    else
    {
      // cv[i] should be split.
      Segment_2    cvi_1, cvi_2;

      seg_traits.curve_split (cv[i],
			      cvi_1, cvi_2,
			      p);

      // The first part should go into c1 and the second into c2.
      c1._push_back (cvi_1);
      c2._push_back (cvi_2);
    }

    // Push all segments labeled (i+1, i+2, ... , n-1) into cv1.
    int    n = cv._size();

    for (j = i+1; j < n; j++)
      c2._push_back (cv[j]);

    return;
  }

  /*!
   * Find the nearest intersection point (or points) of two given curves to
   * the right lexicographically of a given point not including the point
   * itself:
   *  - If the intersection of the two curves is a point to the right of the 
   *    given point, it is returned through both point references p1 and p2.
   *  - If the intersection of the two curves is an X_monotone_curve_2, 
   *    that is, they overlap at infinitely many points, then the rightmost
   *    segment of the intersection is returned.
   * NOTE: When there is an overlap we will always return a SEGMENT (i.e.,
   *       p1 and p2 will be on a segment) even if the overlap is a polyline,
   *       but this is still sufficient for the arrangement.
   * \param cv1 The first curve.
   * \param cv2 The second curve.
   * \param p The refernece point.
   * \param p1 The first output point.
   * \param p2 The second output point.
   * \return (true) if c1 and c2 do intersect to the right of p, or (false)
   * if no such intersection exists.
   */
  bool nearest_intersection_to_right (const X_monotone_curve_2& cv1,
				      const X_monotone_curve_2& cv2,
				      const Point_2& p,
				      Point_2& p1, Point_2& p2) const
  {
    return (_nearest_intersection_to_side (cv1, cv2,
					   p,
					   true,     // To the right of p.
					   p1, p2));
  }

  /*!
   * Find the nearest intersection point (or points) of two given curves to
   * the left lexicographically of a given point not including the point
   * itself:
   *  - If the intersection of the two curves is a point to the left of the 
   *    given point, it is returned through both point references p1 and p2.
   *  - If the intersection of the two curves is an X_monotone_curve_2, 
   *    that is, they overlap at infinitely many points, then the leftmost
   *    segment of the intersection is returned.
   * NOTE: When there is an overlap we will always return a SEGMENT (i.e.,
   *       p1 and p2 will be on a segment) even if the overlap is a polyline,
   *       but this is still sufficient for the arrangement.
   * \param cv1 The first curve.
   * \param cv2 The second curve.
   * \param p The refernece point.
   * \param p1 The first output point.
   * \param p2 The second output point.
   * \return (true) if c1 and c2 do intersect to the right of p, or (false)
   * if no such intersection exists.
   */
  bool nearest_intersection_to_left (const X_monotone_curve_2& cv1,
				     const X_monotone_curve_2& cv2,
				     const Point_2& p,
				     Point_2& p1, Point_2& p2) const
  {
    return (_nearest_intersection_to_side (cv1, cv2,
					   p,
					   false,     // To the left of p.
					   p1, p2));
  }

  /*!
   * Check whether the two given curves overlap.
   * \patam cv1 The first curve.
   * \patam cv2 The second curve.
   * \return (true) if the two curves overlap in a one-dimensional subcurve
   * (i.e., not in a finite number of points). Otherwise, if they have a finite
   * number of intersection points, or none at all, return (false).
   */
  bool curves_overlap(const X_monotone_curve_2& cv1,
                      const X_monotone_curve_2& cv2) const
  {
    // Get the leftmost point of cv1 and of cv2.
    const Point_2& left1 = _is_curve_to_right(cv1) ? curve_source(cv1) : 
                                                     curve_target(cv1);
    const Point_2& left2 = _is_curve_to_right(cv2) ? curve_source(cv2) : 
                                                     curve_target(cv2);

    // Pick the righmost point of the two to start with.
    Point_2 p = (seg_traits.compare_x(left1, left2) == LARGER) ? left1 : left2;

    // Try to find an intersection to the right of p.
    Point_2 p1, p2;

    while (_nearest_intersection_to_side (cv1, cv2,
					  p,
					  true,     // To the right of p.
					  p1, p2))
    {
      // Check if an overlap has been detected:
      if (! seg_traits.point_equal (p1, p2))
	return (true);

      // Otherwise, p1==p2 is the next intersection point to the right.
      // Continue checking the intersections to its right.
      p = p1;
    }

    // If we reached here, no overlap has been detected.
    return (false);
  }

private:

  /*!
   * Return the index of the segment in the polyline cv that contains the
   * point q in its x-range. The function performs a binary search, so if the
   * point q is in the x-range of the polyline with n segments, the segment
   * containing it can be located in O(log n) operations.
   * \param cv The polyline curve.
   * \param q The point.
   * \return An index i such that q is in the x-range of cv[i].
   *         If q is not in the x-range of cv, returns (-1).
   */
  int _locate_point (const X_monotone_curve_2& cv, const Point_2& q) const
  {
    // First check whether the polyline curve really contains q in its x-range.
    int               from = 0;
    Comparison_result res_from;
    int               to = cv._size() - 1;
    Comparison_result res_to;

    res_from = seg_traits.compare_x (seg_traits.curve_source(cv[from]), q);
    if (res_from == EQUAL)
      return (from);
    
    res_to = seg_traits.compare_x (seg_traits.curve_target(cv[to]), q);
    if (res_to == EQUAL)
      return (to);
    
    if (res_from == res_to)
      return (-1);

    // Perform a binary search and locate the segment that contains q in its
    // x-range.
    int               mid;
    Comparison_result res_mid_s, res_mid_t;

    while (to > from)
    {
      mid = (from + to) / 2;

      if (mid > from)
      {
	res_mid_s = seg_traits.compare_x (seg_traits.curve_source(cv[mid]), q);

	if (res_mid_s == EQUAL)
	  return (mid);
	
	if (res_mid_s == res_from)
	  from = mid;
	else
	  to = mid - 1;
      }
      else
      {
	CGAL_assertion (mid < to);
	res_mid_t = seg_traits.compare_x (seg_traits.curve_target(cv[mid]), q);

	if (res_mid_t == EQUAL)
	  return (mid);
	
	if (res_mid_t == res_to)
	  to = mid;
	else
	  from = mid + 1;
      }
    }

    // In case (from == to), and we know that the polyline contains the q:
    CGAL_assertion (from == to);
    return (from);
  }

  /*!
   * Find the index of the segment in the polyline that is defined to the
   * left (or to the right) of a the point q.
   * \param cv The polyline curve.
   * \param q The point.
   * \param to_right (true) if we wish to locate a segment to the right of q,
   *                 (false) if we wish to locate a segment to its right.
   * \return An index i such that cv[i] is defined to the left (or to the
   *         right) of q, or (-1) if no such segment exists.
   */
  int _locate_point_side (const X_monotone_curve_2& cv, const Point_2& q,
			  const bool& to_right) const
  {
    // First locate a segment cv[i] that contains q in its x-range.
    int    i = _locate_point (cv, q);

    if (i < 0)
      return (-1);

    // If we seek an end-point to the right of q, q must be smaller than it.
    // If we seek an end-point to its left, q must be larger.
    const Comparison_result  cres = (to_right) ? SMALLER : LARGER;

    // Check whether x(q) is the x coordinate of an end-point of cv.
    if (seg_traits.compare_x (q, seg_traits.curve_source(cv[i])) == EQUAL)
    {
      // If x(q) is at cv[i]'s source and its target is to the right (or left) 
      // of q, then cv[i] is defined to the right (or left) of q. 
      if (seg_traits.compare_x (q, 
				seg_traits.curve_target(cv[i])) == cres)
      {
	return (i);
      }

      // Otherwise, since cv[i]'s source is cv[i-1]'s target, we should check
      // if the source of the previous curve is to the right (or left) of q.
      if (i > 0 &&
	  seg_traits.compare_x (q, 
			        seg_traits.curve_source(cv[i-1])) == cres)
      {
	return (i-1);
      }
      else
      {
	return (-1);
      }
    }
    else if (seg_traits.compare_x (q, seg_traits.curve_target(cv[i])) == EQUAL)
    {
      // If x(q) is at cv[i]'s target and its source is to the right (or left) 
      // of q, then cv[i] is defined to the right (or left) of q. 
      if (seg_traits.compare_x (q, 
				seg_traits.curve_source(cv[i])) == cres)
      {
	return (i);
      }

      // Otherwise, since cv[i]'s target is cv[i+1]'s source, we should check
      // if the target of the next curve is to the right (or left) of q.
      if (i < (cv._size() - 1) &&
	  seg_traits.compare_x (q, 
			        seg_traits.curve_target(cv[i+1])) == cres)
      {
	return (i+1);
      }
      else
      {
	return (-1);
      }
    }

    // In case q is in cv[i]'s interior:
    return (i);
  }

  /*!
   * Check whether the cuve is defined from the left to the right.
   * \param cv The polyline curve.
   * \return (true) if the curve target is lexicographically larger than its
   *         source; (false) otherwise.
   */
  bool _is_curve_to_right (const X_monotone_curve_2& cv) const
  {
    const Point_2&    ps = seg_traits.curve_source(cv[0]);
    const Point_2&    pt = seg_traits.curve_target(cv[0]);
    Comparison_result res;
    
    res = seg_traits.compare_x (ps, pt);
    if (res == EQUAL)
      res = seg_traits.compare_xy (ps, pt);
				  
    return (res == SMALLER);
  }

  /*!
   * Find the nearest intersection point (or points) of two given curves to
   * the right (or to the left) lexicographically of a given point not 
   * including the point itself.
   * This function is used by both nearest_intersection_to_right() and
   * nearest_intersection_to_left() to avoid code duplication.
   * \param cv1 The first curve.
   * \param cv2 The second curve.
   * \param p The refernece point.
   * \param to_right Should we search for an intersection to the right or to
   *                 the left of p.
   * \param p1 The first output point.
   * \param p2 The second output point.
   * \return (true) if c1 and c2 do intersect to the right of p, or (false)
   * if no such intersection exists.
   */
  bool _nearest_intersection_to_side (const X_monotone_curve_2& cv1,
				      const X_monotone_curve_2& cv2,
				      const Point_2& p,
				      const bool& to_right,
				      Point_2& p1, Point_2& p2) const
  {
    // Get the indices of the segments in cv1 and cv2 containing p.
    int    i1 = _locate_point (cv1, p);
    int    i2 = _locate_point (cv2, p);


    // Check if cv1 and cv2 are defined from left to right, and also
    // determine the desired comparison result (and its inverse).
    Comparison_result d_res;
    Comparison_result i_res;
    int    n1 = cv1._size();
    int    n2 = cv2._size();
    int    inc1;
    int    inc2;

    if (to_right)
    {
      d_res = SMALLER;
      i_res = LARGER;
      inc1 = _is_curve_to_right(cv1) ? 1 : -1;
      inc2 = _is_curve_to_right(cv2) ? 1 : -1;
    }
    else
    {
      d_res = LARGER;
      i_res = SMALLER;
      inc1 = _is_curve_to_right(cv1) ? -1 : 1;
      inc2 = _is_curve_to_right(cv2) ? -1 : 1;
    }

    // Check if the entire curve cv1 is defined to the right (left) of p. 
    if (i1 < 0)
    {
      if (inc1 == 1)
      {
	// The source is the leftmost (rightmost) point in cv1, so p must be 
	// to its left (right).
	// Otherwise, p is to the right (left) of the x-range of cv1. 
	if (seg_traits.compare_x(p, 
				 seg_traits.curve_source(cv1[0])) == d_res)
	  i1 = 0;
	else
	  return (false);
      }
      else // if (inc1 == -1)
      {
	// The target is the leftmost (rightmost) point in cv1, so p must be 
	// to its left (right).
	// Otherwise, p is to the right (left) of the x-range of cv1. 
	if (seg_traits.compare_x(p,
				 seg_traits.curve_target(cv1[n1-1])) == d_res)
	  i1 = n1-1;
	else
	  return (false);
      }
    }    

    // Check if the entire curve cv2 is defined to the right (left) of p. 
    if (i2 < 0)
    {
      if (inc2 == 1)
      {
	// The source is the leftmost (rightmost) point in cv2, so p must be 
	// to its left (right).
	// Otherwise, p is to the right (left) of the x-range of cv2. 
	if (seg_traits.compare_x(p, 
				 seg_traits.curve_source(cv2[0])) == d_res)
	  i2 = 0;
	else
	  return (false);
      }
      else // if (inc2 == -1)
      {
	// The target is the leftmost (rightmost) point in cv2, so p must be 
	// to its left (right).
	// Otherwise, p is to the right (left) of the x-range of cv2. 
	if (seg_traits.compare_x(p,
				 seg_traits.curve_target(cv2[n2-1])) == d_res)
	  i2 = n2-1;
	else
	  return (false);
      }
    }

    // Try to locate the intersection point.
    bool               found;
    Comparison_result  res;

    while (i1 >= 0 && i1 < n1 && i2 >= 0 && i2 < n2)
    {
      // Check if the two current segment intersect to the right (left) of p.
      if (to_right)
	found = seg_traits.nearest_intersection_to_right (cv1[i1], cv2[i2],
							  p,
							  p1, p2);
      else
	found = seg_traits.nearest_intersection_to_left (cv1[i1], cv2[i2],
							  p,
							  p1, p2);

      if (found)
      {
	// In case an overlap was detected, stop here:
	if (! seg_traits.point_equal (p1, p2))
	  return (true);

	// In case we found a single intersection point, check whether it
	// is the next end-point of cv1[i1] or of cv2[i2].
	bool    eq1, eq2;

	eq1 = seg_traits.point_equal (p1, 
				      (inc1 == 1) ? 
				      seg_traits.curve_target(cv1[i1]) :
				      seg_traits.curve_source(cv1[i1]));

	eq2 = seg_traits.point_equal (p1,
				      (inc2 == 1) ? 
				      seg_traits.curve_target(cv2[i2]) :
				      seg_traits.curve_source(cv2[i2]));

	if (!eq1 && !eq2)
	  return (true);

	// Proceed to the next curves.
	if (eq1)
	  i1 += inc1;
	if (eq2)
	  i2 += inc2;

	if (i1 < 0 || i1 >= n1 || i2 < 0 || i2 >= n2)
	  return (true);

	// Check if the next curves overlap, and the nearest overlap endpoint
	// equals p1.
	Point_2  q1, q2;

	if (to_right)
	  found = seg_traits.nearest_intersection_to_right (cv1[i1], cv2[i2],
							    p,
							    q1, q2);
	else
	  found = seg_traits.nearest_intersection_to_left (cv1[i1], cv2[i2],
							   p,
							   q1, q2);


	if (found &&
	    ! seg_traits.point_equal (q1, q2) &&
	    seg_traits.point_equal (p1, q1))
	{
	  // Now p1 (== q1) --> q2 is an overlapping segment.
	  p2 = q2;
	}

	return (true);
      }

      // Find the segment whose end-point is the leftmost (rightmost) and move
      //  forward (or backward) on its polyline.
      res = seg_traits.compare_x ((inc1 == 1) ? 
				  seg_traits.curve_target(cv1[i1]) :
				  seg_traits.curve_source(cv1[i1]),
				  (inc2 == 1) ? 
				  seg_traits.curve_target(cv2[i2]) :
				  seg_traits.curve_source(cv2[i2]));

      if (res == d_res)
      {
	i1 += inc1;
      }
      else if (res == i_res)
      {
	i2 += inc2;
      }
      else
      {
	i1 += inc1;
	i2 += inc2;
      }
    }

    // No intersection found:
    return (false);
  }

};

/*!
 * A representation of a polyline, as used by the Arr_polyline_traits_2
 * traits class.
 */
template <class Segment_traits_>
class Polyline_2
{
  friend class Arr_polyline_traits_2<Segment_traits_>;

public:

  typedef Segment_traits_                           Segment_traits_2;
  typedef typename Segment_traits_2::Point_2        Point_2;
  typedef typename Segment_traits_2::Curve_2        Segment_2;
  typedef ::std::vector<Segment_2>                  Base;

private:
  
  // The segments that comprise the poyline:
  ::std::vector<Segment_2>                          segments;

public:

  /*!
   * Default constructor.
   */
  Polyline_2 () :
    segments()
  {}

  /*!
   * Constructor from a range of points, defining the endpoints of the
   * polyline segments.
   * \param pts_begin An iterator pointing to the first point in the range.
   * \param pts_end An iterator pointing after the last point in the range.
   * \pre The are at least 2 points in the range.
   *      In other cases, an empty polyline will be created.
   */
  template <class InputIterator>
  Polyline_2 (const InputIterator& pts_begin,
	      const InputIterator& pts_end) :
    segments()
  {
    // Check if there are no points in the range:
    InputIterator  ps = pts_begin;

    if (ps == pts_end)
      return;

    // Construct a segment from each to adjacent points.
    InputIterator  pt = ps;
    pt++;

    while (pt != pts_end)
    {
      segments.push_back (Segment_2 (*ps, *pt));
      ps++; pt++;
    }
  }

  /*!
   * Create a bounding-box for the polyline.
   * \return The bounding-box.
   */
  Bbox_2 bbox() const
  {
    // Compute the union of the bounding boxes of all segments.
    Bbox_2 bbox;
    int    n = _size();
    int    i;

    for (i = 0; i < n; i++)
    {
      if (i == 0)
	bbox = (*this)[i].bbox();
      else
	bbox = bbox + (*this)[i].bbox();
    }

    return (bbox);
  }

  class const_iterator;
  friend class const_iterator;

  /*!
   * An iterator for the polyline points.
   */
  class const_iterator
  {
  private:

    const Polyline_2<Segment_traits_> *cvP;        // The polyline curve.
    int                               n_pts;       // Its number of points.
    Segment_traits_                   seg_traits;  // Auxiliary variable.
    int                               i;           // The current point.

    /*!
     * Private constructor.
     * \param cv The scanned curve.
     * \param index The index of the segment.
     */
    const_iterator (const Polyline_2<Segment_traits_>* _cvP,
		    const int& index) :
      cvP(_cvP),
      i(index)
    {
      if (cvP == NULL)
	n_pts = 0;
      else
	n_pts = (cvP->_size() == 0) ? 0 : (cvP->_size() + 1);
    }

  public:
    
    /*!
     * Default constructor.
     */
    const_iterator () :
      cvP(NULL),
      n_pts(0),
      i(-1)
    {}

    /*!
     * Dereference operator.
     * \return The current point.
     */
    Point_2 operator* () const
    {
      CGAL_assertion(cvP != NULL);
      CGAL_assertion(i >= 0 && i < n_pts);

      if (i == 0)
	// First point is the source of the first segment.
	return (seg_traits.curve_source ((*cvP)[0]));
      else
	// Return the target of the (i-1)'st segment
	return (seg_traits.curve_target ((*cvP)[i-1]));
    }

    /*!
     * Increment operators.
     */
    void operator++ () 
    {
      if (cvP != NULL && i < n_pts)
	i++;
    }

    void operator++ (int)
    {
      if (cvP != NULL && i < n_pts)
	i++;
    }

    /*!
     * Decrement operators.
     */
    void operator-- ()
    {
      if (cvP != NULL && i >= 0)
	i--;
    }

    void operator-- (int)
    {
      if (cvP != NULL && i >= 0)
	i--;
    }

    /*!
     * Equality operators.
     */
    bool operator== (const const_iterator& it) const
    {
      return (cvP == it.cvP && i == it.i);
    }

    bool operator!= (const const_iterator& it) const
    {
      return (cvP != it.cvP || i != it.i);
    }

    friend class Polyline_2<Segment_traits_>;
  };

  /*!
   * Get an iterator for the polyline points.
   * \return An iterator that points on the first point.
   */
  const_iterator begin () const
  {
    if (_size() == 0)
      return (const_iterator (NULL, -1));
    else
      return (const_iterator (this, 0));
  }

  /*!
   * Get a past-the-end iterator for the polyline points.
   * \return A past-the-end iterator.
   */
  const_iterator end () const
  {
    if (_size() == 0)
      return (const_iterator (NULL, -1));
    else
      return (const_iterator (this, _size() + 1));
  }

  /*!
   * Get an reverse iterator for the polyline points.
   * \return An iterator that points on the last point.
   */
  const_iterator rbegin () const
  {
    if (_size() == 0)
      return (const_iterator (NULL, -1));
    else
      return (const_iterator (this, _size()));
  }

  /*!
   * Get a reverse past-the-end iterator for the polyline points.
   * \return A reverse past-the-end iterator.
   */
  const_iterator rend () const
  {
    if (_size() == 0)
      return (const_iterator (NULL, -1));
    else
      return (const_iterator (this, -1));
  }

  /*!
   * Get the number of points contained in the polyline.
   * \return The number of points.
   */
  unsigned int points () const
  {
    if (_size() == 0)
      return (0);
    else
      return (_size() + 1);
  }

private:

  /*!
   * Append a segment to the polyline.
   * \param seg The new segment to be appended to the polyline.
   * \pre If the polyline is not empty, the segment source must be the
   *      same as the target point of the last segment in the polyline.
   */
  inline void _push_back (const Segment_2& seg)
  {
    segments.push_back(seg);
    return;
  }

  /*!
   * Get the number of segments that comprise the poyline.
   * \return The number of segments.
   */
  inline int _size () const
  {
    return (static_cast<int>(segments.size()));
  }

  /*!
   * Get the i'th segment of the polyline.
   * \param i The segment index (from 0 to size()-1).
   * \return A const reference to the segment.
   */
  inline const Segment_2& operator[] (const int& i) const
  {
    return (segments[i]);
  }

  /*!
   * Clear the polyline.
   */
  inline void _clear ()
  {
    segments.clear();
    return;
  }
};

/*!
 * Output operator for a polyline.
 */
template <class Segment_traits_, class Stream_>
Stream_ & operator<< (Stream_ & os, const Polyline_2<Segment_traits_> & cv)
{
  typename Polyline_2<Segment_traits_>::const_iterator ps = cv.begin();
  typename Polyline_2<Segment_traits_>::const_iterator pt = ps; pt++;

  while (pt != cv.end()) {
    typename Segment_traits_::Curve_2 seg(*ps, *pt);
    os << seg;
    ps++; pt++;
  }
  return (os);
}

CGAL_END_NAMESPACE

#endif
