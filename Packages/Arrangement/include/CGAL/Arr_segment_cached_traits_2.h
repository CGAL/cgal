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
// file          : include/CGAL/Pm_segment_traits_2.h
// package       : Planar_map (5.87)
// maintainer    : Eyal Flato        <flato@math.tau.ac.il>
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Efi Fogel         <efif@post.tau.ac.il>
//                 Ron Wein          <wein@post.tau.ac.il>
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ======================================================================
#ifndef CGAL_ARR_SEGMENT_CACHED_TRAITS_2_H
#define CGAL_ARR_SEGMENT_CACHED_TRAITS_2_H

#include <CGAL/Planar_map_2/Pm_segment_utilities_2.h>
#include <CGAL/tags.h>
#include <CGAL/intersections.h>

#include <list>
#include <fstream>

CGAL_BEGIN_NAMESPACE

template <class Kernel_>
class Arr_segment_cached_traits_2 : public Kernel_
{
public:
  typedef Kernel_                         Kernel;

  typedef Tag_true                        Has_left_category;

  /*!
   * Representation of a segement with cached data.
   */
  class Segment_cached_2
  {
    typedef typename Kernel_::Segment_2             Segment_2;
    typedef typename Kernel_::Point_2               Point_2;

  private:

    bool      is_orig;          // Is this an original segment.
    Segment_2 orig_seg;         // The original segment.
    Point_2   ps, pt;           // The source a target points.
    bool      is_vert;          // Is this a vertical segment.

  public:

    /*!
     * Default constructor.
     */
    Segment_cached_2 () :
      is_orig(true),
      orig_seg(),
      is_vert(-1)
    {}

    /*!
     * Constructor from a segment.
     */
    Segment_cached_2 (const Segment_2& seg) :
      is_orig(true),
      orig_seg(seg)
    {
      Kernel_   kernel;
      
      is_vert = kernel.is_vertical_2_object()(seg);
      
      typename Kernel_::Construct_vertex_2 
	construct_vertex = kernel.construct_vertex_2_object();

      ps = construct_vertex(seg, 0);
      pt = construct_vertex(seg, 1);
    }

    /*!
     * Constructor from two points.
     */
    Segment_cached_2 (const Point_2& source, const Point_2& target) :
      is_orig(true),
      orig_seg(source,target)
    {
      Kernel_   kernel;
      
      is_vert = kernel.is_vertical_2_object()(orig_seg);
      
      ps = source;
      pt = target;
    }

    /*!
     * Cast to a segment.
     */
    operator Segment_2 () const
    {
      return (Segment_2(ps, pt));
    }

    /*!
     *
     */
    CGAL::Bbox_2 bbox()
    {
      Segment_2 seg(ps, pt);
      return seg.bbox();
    }
      
  private:

    friend class Arr_segment_cached_traits_2;
  };

  // Traits objects
  typedef typename Kernel::Point_2        Point_2;
  typedef Segment_cached_2                X_curve_2;
  typedef Segment_cached_2                Curve_2;

  // Backward compatability    
  typedef Point_2                         Point;
  typedef X_curve_2                       X_curve;
  typedef X_curve_2                       Curve;

  // Currently, I leave this in the traits
  // Maybe we can change the usage inside Planar_map_2
  typedef enum
  {
    UNDER_CURVE        = -1,
    CURVE_NOT_IN_RANGE =  0,
    ABOVE_CURVE        =  1,
    ON_CURVE           =  2

  } Curve_point_status;	

protected:
  // Functors:
  typedef typename Kernel::Is_vertical_2        Is_vertical_2;
  typedef typename Kernel::Construct_vertex_2   Construct_vertex_2;
  typedef typename Kernel::Less_x_2             Less_x_2;
  typedef typename Kernel::Equal_2              Equal_2;
  typedef typename Kernel::Compare_x_2          Compare_x_2;
  typedef typename Kernel::Compare_xy_2         Compare_xy_2;
    
  typedef CGAL::Counterclockwise_in_between_for_segments_2<Kernel, X_curve_2>
                                                Counterclockwise_in_between_2;

protected:
  inline Counterclockwise_in_between_2 counterclockwise_in_between_2_object()
    const
  { return Counterclockwise_in_between_2(); }

public:
  // Creation
  Arr_segment_cached_traits_2() {}

  // Operations
  // ----------
    
  /*! compare_x() compares the x-coordinates of two given points
   * \param p1 the first point
   * \param p2 the second point
   * \return LARGER if x(p1) > x(p2), SMALLER if x(p1) < x(p2), or else EQUAL
   *
   * \todo replace indirect use compare_x() with compare_x_2()
   */
  Comparison_result compare_x(const Point_2 & p1, const Point_2 & p2) const
  {
    return compare_x_2_object()(p1, p2);
  }

  /*! compare_xy() compares lexigoraphically the two points by x, then by y.
   * \param p1 the first point
   * \param p2 the second point
   * \return LARGER if x(p1) > x(p2), or if x(p1) = x(p2) and y(p1) > y(p2); 
   *         SMALLER if x(p1) < x(p2), or if x(p1) = x(p2) and y(p1) < y(p2);
   *         or else EQUAL
   */
  Comparison_result compare_xy(const Point_2 & p1, const Point_2 & p2) const
  {
    return compare_xy_2_object()(p1, p2);
  }

  /*! curve_is_vertical()
   * \param cv the curve
   * \return true iff the curve is vertical
   */
  bool curve_is_vertical(const X_curve_2 & cv) const 
  {
    return (cv.is_vert);
  } 

  /*! curve_is_in_x_range()
   * \param cv the curve
   * \param q the point
   * \return true if q is in the x range of cv
   */
  bool curve_is_in_x_range(const X_curve_2 & cv, const Point_2 & q) const
  {
    Compare_x_2       compare_x = compare_x_2_object();
    Comparison_result res1 = compare_x (q, cv.ps);
    Comparison_result res2 = compare_x (q, cv.pt);

    // We check if x(p) equals the x value of one of the end-points.
    // If not, we check whether one end-point is to p's left and the other is
    // to its right.
    return ((res1 == EQUAL) || (res2 == EQUAL) ||
	    (res1 != res2));
  }

  /*! curve_compare_at_x() compares the y-coordinate of two given curves at
   * the x-coordinate of a given point.
   * Preconditions: The point q is in the x range of the two curves.
   * \param cv1 the first curve
   * \param cv2 the second curve
   * \param q the point
   * \return LARGER if cv1(x(q)) > cv2(x(q)), SMALLER if cv1(x(q)) < cv2(x(q),
   *  or else EQUAL.
   *
   * \todo replace indirect use curve_compare_at_x() with compare_y_at_x_2()
   */
  Comparison_result curve_compare_at_x(const X_curve_2 & cv1, 
				       const X_curve_2 & cv2, 
				       const Point_2 & q) const
  {
    CGAL_precondition(curve_is_in_x_range(cv1, q));
    CGAL_precondition(curve_is_in_x_range(cv2, q));

    // Compare using the original segments.
    return (compare_y_at_x_2_object()(q, cv1.orig_seg, cv2.orig_seg));
  }

  /*! curve_compare_at_x_left() compares the y value of two curves in an
   * epsilon environment to the left of the x value of the input point
   * Preconditions: The point q is in the x range of the two curves, and both
   * of them must be also be defined to its left.
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
    
    // Since the curves are continuous, if they are not equal at q, the same
    // result also applies to q's left.
    Comparison_result res = 
      compare_y_at_x_2_object()(q, cv1.orig_seg, cv2.orig_seg);

    if (res != EQUAL)
      return (res);     
    
    // <cv2> and <cv1> meet at a point with the same x-coordinate as q
    // compare their derivatives.
    // Notice we use the original segments in order to compare the slopes.
    return (compare_slope_2_object()(cv2.orig_seg, cv1.orig_seg));
  }
    
  /*! curve_compare_at_x_right() compares the y value of two curves in an
   * epsilon environment to the right of the x value of the input point
   * Preconditions: The point q is in the x range of the two curves, and both
   * of them must be also be defined to its right.
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
    
    // Since the curves are continuous, if they are not equal at q, the same
    // result also applies to q's right.
    Comparison_result res =
      compare_y_at_x_2_object()(q, cv1.orig_seg, cv2.orig_seg);

    if (res != EQUAL)
      return (res);     
    
    // <cv1> and <cv2> meet at a point with the same x-coordinate as q
    // compare their derivatives
    // Notice we use the original segments in order to compare the slopes.
    return (compare_slope_2_object()(cv1.orig_seg, cv2.orig_seg));
  }
    
  /*! 
   * Return the curve-point status of the input objects
   */
  Curve_point_status curve_get_point_status (const X_curve_2 & cv, 
					     const Point_2 & p) const
  {
    if (! curve_is_in_x_range(cv, p))
      return (CURVE_NOT_IN_RANGE);

    if (! cv.is_vert)
    {
      // Compare with the original segment.
      Comparison_result res = compare_y_at_x_2_object()(p, cv.orig_seg);
      return ((res == LARGER) ? ABOVE_CURVE :
	      ((res == SMALLER) ? UNDER_CURVE : ON_CURVE));
    }
    else
    {
      // Compare with the vertical segment's end-points.
      Compare_xy_2      compare_xy = compare_xy_2_object();
      Comparison_result res1 = compare_xy (p, cv.ps);
      Comparison_result res2 = compare_xy (p, cv.pt);
      
      if (res1 == LARGER && res2 == LARGER)
	return (ABOVE_CURVE);
      else if (res1 == SMALLER && res2 == SMALLER)
	return (UNDER_CURVE);
      else
	return (ON_CURVE);
    }
  }

  /*! \todo Degenerate cases may not work! Talk with Eyal to fix the actual
   * code in Pmwx to use the same consisting definitions of
   * curve_is_between_cw(), counterclockwise_in_between_2_object(), and
   * the kernel function that is used to implement the later.
   */
  bool curve_is_between_cw(const X_curve_2 & cv, 
                           const X_curve_2 & first, 
                           const X_curve_2 & second, 
                           const Point_2 & point) const
  {
    // Since I'm lazy:
    typename Kernel::Segment_2   cvx (cv.ps, cv.pt);
    typename Kernel::Segment_2   cv1 (first.ps, first.pt);
    typename Kernel::Segment_2   cv2 (second.ps, second.pt);

    return counterclockwise_in_between_2_object()(point, cvx, cv1, cv2);
    // Notice the change in order of first and second
    //return counterclockwise_in_between_2_object()(point, cv, second, first);
  }

  /*! 
   * Check if the two curves are the same (have the same graph).
   */
  bool curve_is_same(const X_curve_2 & cv1,const X_curve_2 & cv2) const
  {
    Equal_2 equal = equal_2_object();
    return ((equal(cv1.ps, cv2.ps) && equal(cv1.pt, cv2.pt)) ||
	    (equal(cv1.ps, cv2.pt) && equal(cv1.pt, cv2.ps)));
  }

  /*!
   * Check if the two points are the same.
   */
  bool point_is_same(const Point_2 & p1, const Point_2 & p2) const
  {
    return equal_2_object()(p1, p2);
  }
  
  /*!
   * Get the curve source.
   */
  const Point_2& curve_source(const X_curve_2 & cv) const 
  { 
    return (cv.ps);
  }

  /*!
   * Get the curve target.
   */
  const Point_2& curve_target(const X_curve_2 & cv) const 
  { 
    return (cv.pt);
  }

  /*!
   * Return TRUE, since a sgement is always x-monotone.
   */
  bool is_x_monotone(const Curve_2 &) const
  {
    return true;
  }
  
  /*!
   * Do nothing, since a sgement is always x-monotone.
   */
  void make_x_monotone(const Curve_2&, std::list<Curve_2>& ) const
  {} 

  /*!
   * Flip a given curve.
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
   * \pre split_pt is on cv but is not an endpoint.
   */
  void curve_split(const X_curve_2& cv, 
		   X_curve_2& c1, X_curve_2& c2, 
                   const Point_2& p) const
  {
    // Check preconditions.
    CGAL_precondition(curve_get_point_status(cv, p) == ON_CURVE);
    CGAL_precondition_code(Compare_xy_2 compare_xy = compare_xy_2_object());
    CGAL_precondition(compare_xy(cv.ps, p) != EQUAL);
    CGAL_precondition(compare_xy(cv.pt, p) != EQUAL);
    
    // Do the split.
    c1.is_orig = false;
    c1.orig_seg = cv.orig_seg;
    c1.ps = cv.ps;
    c1.pt = p;
    c1.is_vert = cv.is_vert;

    c2.is_orig = false;
    c2.orig_seg = cv.orig_seg;
    c2.ps = p;
    c2.pt = cv.pt;
    c2.is_vert = cv.is_vert;

    return;
  }

  /*! 
   * Compares the location of the intersection point of two given curves 
   * with a given point.
   * \param cv1 the first curve
   * \param cv2 the second curve
   * \param p the point
   * \return true if cv1 and cv2 intersect at a point that is lexicographically
   * larger than p, that is above or to the right of p but not on p.    
   */
  //returns true iff the intersection is strictly right of pt
  bool do_intersect_to_right(const X_curve_2& cv1, const X_curve_2& cv2,
                             const Point_2& p) const 
  {
    bool     is_overlap;
    Point_2  ip1, ip2;

    if (! _find_intersection (cv1, cv2, is_overlap, ip1, ip2))
      return (false);

    if (! is_overlap) 
    {
      return (compare_xy_2_object()(ip1, p) == LARGER);
    }
    else
    {
      // Since always ip1 < ip2.
      return (compare_xy_2_object()(ip2, p) == LARGER);
    }
  }

  /*! do_intersect_to_left() compares the location of the intersection
   * point of two given curves with a given point.
   * \param cv1 the first curve
   * \param cv2 the second curve
   * \param p the point
   * \return true if c1 and c2 intersect at a point that is lexicographically
   * smaller than pt, that is bellow or to the left of pt but not on pt.    
   */
  bool do_intersect_to_left(const X_curve_2 & cv1, const X_curve_2 & cv2,
                            const Point_2 & p) const 
  {
    bool     is_overlap;
    Point_2  ip1, ip2;

    if (! _find_intersection (cv1, cv2, is_overlap, ip1, ip2))
      return (false);

    if (! is_overlap) 
    {
      return (compare_xy_2_object()(ip1, p) == SMALLER);
    }
    else
    {
      // Since always ip1 < ip2.
      return (compare_xy_2_object()(ip1, p) == SMALLER);
    }
  }

  /*! nearest_intersection_to_right() finds the nearest intersection point of
   * two given curves to the right of a given point. Nearest is defined as the
   * lexicographically nearest not including the point itself with one
   * exception explained bellow..
   * If the intersection of the two curves is an X_curve_2, that is,
   * there is an overlapping subcurve, then if the the source and target of the
   * subcurve are strickly to the right, they are returned through two
   * other point references p1 and p2. If pt is between the source and target
   * of the overlapping subcurve, or pt is its left endpoint, pt and the target
   * of the right endpoint of the subcurve are returned through p1 and p2 
   * respectively.
   * If the intersection of the two curves is a point to the right of pt, pt
   * is returned through the p1 and p2.
   * \param cv1 the first curve
   * \param cv2 the second curve
   * \param p the point to compare against
   * \param p1 the first point reference
   * \param p2 the second point reference
   * \return true if c1 and c2 do intersect to the right of pt. Otherwise,
   * false
   */
  bool nearest_intersection_to_right(const X_curve_2 & cv1,
                                     const X_curve_2 & cv2,
                                     const Point_2 & p,
                                     Point_2 & p1, Point_2 & p2) const
  {
    bool     is_overlap;
    Point_2  ip1, ip2;

    if (! _find_intersection (cv1, cv2, is_overlap, ip1, ip2))
      return (false);

    if (! is_overlap) 
    {
      if (compare_xy_2_object()(ip1, p) == LARGER)
      {
	p1 = ip1;
	p2 = ip1;
	return (true);
      }
      return (false);
    }
    else
    {
      // Notice that ip1 < ip2.
      if (compare_xy_2_object()(ip1, p) == LARGER)
      {
	p1 = ip1;
	p2 = ip2;
	return (true);
      }
      else if (compare_xy_2_object()(ip2, p) == LARGER)
      {
	p1 = p;
	p2 = ip2;
	return (true);
      }
      return (false);
    }
  }

  /*! nearest_intersection_to_left() finds the nearest intersection point of
   * two given curves to the left of a given point. Nearest is defined as the
   * lexicographically nearest not including the point itself with one
   * exception explained bellow..
   * If the intersection of the two curves is an X_curve_2, that is,
   * there is an overlapping subcurve, then if the the source and target of the
   * subcurve are strickly to the left, they are returned through two
   * other point references p1 and p2. If pt is between the source and target
   * of the overlapping subcurve, or pt is its left endpoint, pt and the target
   * of the left endpoint of the subcurve are returned through p1 and p2 
   * respectively.
   * If the intersection of the two curves is a point to the left of pt, pt
   * is returned through the p1 and p2.
   * \param cv1 the first curve
   * \param cv2 the second curve
   * \param p the point to compare against
   * \param p1 the first point reference
   * \param p2 the second point reference
   * \return true if c1 and c2 do intersect to the left of pt. Otherwise,
   * false
   */
  bool nearest_intersection_to_left(const X_curve_2 & cv1,
                                    const X_curve_2 & cv2,
                                    const Point_2 & p,
                                    Point_2 & p1, Point_2 & p2) const
  {
    bool     is_overlap;
    Point_2  ip1, ip2;

    if (! _find_intersection (cv1, cv2, is_overlap, ip1, ip2))
      return (false);

    if (! is_overlap) 
    {
      if (compare_xy_2_object()(ip1, p) == SMALLER)
      {
	p1 = ip1;
	p2 = ip1;
	return (true);
      }
      return (false);
    }
    else
    {
      // Notice that ip1 < ip2.
      if (compare_xy_2_object()(ip2, p) == SMALLER)
      {
	p1 = ip1;
	p2 = ip2;
	return (true);
      }
      else if (compare_xy_2_object()(ip1, p) == SMALLER)
      {
	p1 = ip1;
	p2 = p;
	return (true);
      }
      return (false);
    }
  }

  /*! curves_overlap() test overlapping between two given curves
   * \patam cv1 the first curve
   * \patam cv2 the second curve
   * \return true if c1 and c2 overlap in a one-dimensional subcurve
   * (i.e., not in a finite number of points). Otherwise, false.
   * \todo end point coincidence instead of intersection!
   */
  bool curves_overlap(const X_curve_2 & cv1, const X_curve_2 & cv2) const
  {
    bool     is_overlap;
    Point_2  ip1, ip2;

    if (! _find_intersection (cv1, cv2, is_overlap, ip1, ip2))
      return (false);

    return (is_overlap); 
  }

private:

  bool _find_intersection (const X_curve_2& cv1,
			   const X_curve_2& cv2,
                           bool& is_overlap,
			   Point_2& p1, Point_2& p2) const
  {
    // Intersect the two original segments.
    Point_2   ip;
    X_curve_2 iseg;
    Object    res = intersect_2_object()(cv1.orig_seg, cv2.orig_seg);

    if (assign(ip, res)) 
    {
      is_overlap = false;

      // Simple case of intersection at a single point.
      if (curve_get_point_status(cv1, ip) == ON_CURVE &&
	  curve_get_point_status(cv2, ip) == ON_CURVE)
      {
	p1 = ip;
	return (true);
      }
      
      return (false);
    }
    if (assign(iseg, res)) 
    {
      // Not implemented yet!!!!
      CGAL_assertion (false);
      // Make sure that at the end, ip1 < ip2 !
    }

    // No intersection at all:
    return (false);
  }

};

template <class Kernel_>
std::ostream& operator<< 
  (std::ostream& os, 	
   const typename Arr_segment_cached_traits_2<Kernel_>::Segment_cached_2& seg)
{
  os << static_cast<typename Kernel_::Segment_2>(seg);
  return (os);
}

CGAL_END_NAMESPACE

#endif
