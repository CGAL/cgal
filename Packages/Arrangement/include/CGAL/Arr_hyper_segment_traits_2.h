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
// file          : include/CGAL/Arr_hyper_segment_traits_2.h
// package       : Arrangement_2 (5.87)
// maintainer    : Efi Fogel         <efif@post.tau.ac.il>
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Ron Wein          <wein@post.tau.ac.il>
//
// coordinator   : Tel-Aviv University (Dan Halperin <danha@post.tau.ac.il>)
//
// ======================================================================
#ifndef CGAL_ARR_HYPER_SEGMENT_TRAITS_2_H
#define CGAL_ARR_HYPER_SEGMENT_TRAITS_2_H

#include <CGAL/tags.h>
#include <CGAL/Arrangement_2/Hyper_segment_2.h>

CGAL_BEGIN_NAMESPACE

/*!
 * A traits class for maintaining an arrangement of x-monotone segments of 
 * canonical hyperbolas.
 */
template <class Kernel_>
class Arr_hyper_segment_traits_2 : public Kernel_
{
public:
  typedef Kernel_                         Kernel;

  typedef Tag_true                        Has_left_category;

  // Traits objects:
  typedef typename Kernel::Point_2        Point_2;
  typedef Hyper_segment_2<Kernel>         X_monotone_curve_2;
  typedef Hyper_segment_2<Kernel>         Curve_2;

  // For backward compatability:    
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

public:

  /*!
   * Default constructor.
   */
  Arr_hyper_segment_traits_2() {}
    
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
    // Currently no vertical segments are allowed:
    return (false);
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
    return (cv.point_is_in_x_range(q));
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

    return (cv1.compare_y_at_x (cv2, q));
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
    // The two curves must not be vertical.
    CGAL_precondition(! curve_is_vertical(cv1));
    CGAL_precondition(! curve_is_vertical(cv2));

    // The two curve must be defined at q and also to its left.
    CGAL_precondition_code(Less_x_2 less_x = less_x_2_object(););

    CGAL_precondition (point_in_x_range(cv1, q));
    CGAL_precondition (less_x(cv1.source(), q) || less_x(cv1.target(), q));
    
    CGAL_precondition (point_in_x_range(cv2, q));
    CGAL_precondition (less_x(cv2.source(), q) || less_x(cv2.target(), q));
    
    // Notice q is a placeholder for the x coordinate of the two curves.
    // That is, if we compare them at x(q) the result should be EQUAL.
    CGAL_precondition(cv1.compare_y_at_x (cv2, q) == EQUAL);
    
    // Compare the slopes of the two segments to determine thir relative
    // position immediately to the left of q.
    // Notice we use the supporting lines in order to compare the slopes.
    return (cv2.compare_slopes (cv1, q));
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
    // The two curves must not be vertical.
    CGAL_precondition(! curve_is_vertical(cv1));
    CGAL_precondition(! curve_is_vertical(cv2));

    // The two curve must be defined at q and also to its right.
    CGAL_precondition_code(Less_x_2 less_x = less_x_2_object(););

    CGAL_precondition (point_in_x_range(cv1, q));
    CGAL_precondition (less_x(q, cv1.source()) || less_x(q, cv1.target()));
    
    CGAL_precondition (point_in_x_range(cv2, q));
    CGAL_precondition (less_x(q, cv2.source()) || less_x(q, cv2.target()));
    
    // Notice q is a placeholder for the x coordinate of the two curves.
    // That is, if we compare them at x(q) the result should be EQUAL.
    CGAL_precondition(cv1.compare_y_at_x (cv2, q) == EQUAL);
    
    // Compare the slopes of the two segments to determine thir relative
    // position immediately to the left of q.
    // Notice we use the supporting lines in order to compare the slopes.
    return (cv1.compare_slopes (cv2, q));
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

    return (cv.point_position(p));
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
    return (cv1.is_equal (cv2));
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
  const Point_2& curve_source(const X_monotone_curve_2 & cv) const 
  { 
    return (cv.source());
  }

  /*!
   * Get the curve target.
   * \param cv The curve.
   * \return The target point.
   */
  const Point_2& curve_target(const X_monotone_curve_2 & cv) const 
  { 
    return (cv.target());
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
    // Return true, since a hyper-segment is always x-monotone.
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
  OutputIterator curve_make_x_monotone (const Curve_2& cv,
					OutputIterator oi) const
  {
    *oi++ = cv;
    return (oi);
  } 

  /*!
   * Flip a given curve.
   * \param cv The input curve.
   * \return The flipped curve. In case of segments, if the input is [s,t],
   * then the flipped curve is simply [t,s].
   */
  X_monotone_curve_2 curve_opposite (const X_monotone_curve_2 & cv) const
  {
    return (cv.flip());
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
  void curve_split(const X_monotone_curve_2& cv, 
		   X_monotone_curve_2& c1, X_monotone_curve_2& c2, 
                   const Point_2& p) const
  {
    cv.split (p, c1, c2);
    return;
  }

  /*!
   * Find the nearest intersection point (or points) of two given curves to
   * the right lexicographically of a given point not includin the point
   * itself, (with one exception explained below).
   * If the intersection of the two curves is an X_monotone_curve_2, that is,
   * they overlap at infinitely many points, then if the right endpoint and the
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
  bool nearest_intersection_to_right(const X_monotone_curve_2 & cv1,
                                     const X_monotone_curve_2 & cv2,
                                     const Point_2 & p,
                                     Point_2 & p1, Point_2 & p2) const
  {
    int      n_pts = cv1.intersect(cv2, p1, p2);

    if (n_pts == 0)
    {
      return (false);
    }
    else if (n_pts == 1) 
    {
      // Check if the intersection is to p's right.
      if (compare_xy_2_object()(p1, p) == LARGER)
      {
	p2 = p1;
	return (true);
      }
      else
	return (false);
    }
    else // In case of an overlap.
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
   * If the intersection of the two curves is an X_monotone_curve_2, that is,
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
  bool nearest_intersection_to_left(const X_monotone_curve_2 & cv1,
                                    const X_monotone_curve_2 & cv2,
                                    const Point_2 & p,
                                    Point_2 & p1, Point_2 & p2) const
  {
    int      n_pts = cv1.intersect(cv2, p1, p2);

    if (n_pts == 0)
    {
      return (false);
    }
    else if (n_pts == 1) 
    {
      // Check if the intersection is to p's left.
      if (compare_xy_2_object()(p1, p) == SMALLER)
      {
	p2 = p1;
	return (true);
      }
      return (false);
    }
    else // In case of an overlap.
    {
      // Notice that p1 < p2.
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
  bool curves_overlap(const X_monotone_curve_2 & cv1,
                      const X_monotone_curve_2 & cv2) const
  {
    Point_2  p1, p2;

    return (cv1.intersect(cv2, p1, p2) == 2);
  }
};

CGAL_END_NAMESPACE

#endif
