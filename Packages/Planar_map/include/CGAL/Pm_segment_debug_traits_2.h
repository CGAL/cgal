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
// Author(s)     : Iddo Hanniel      <hanniel@math.tau.ac.il>
//                 Eyal Flato        <flato@post.tau.ac.il>
//                 Oren Nechushtan   <theoren@math.tau.ac.il>
//                 Eti Ezra          <estere@post.tau.ac.il>
//                 Shai Hirsch       <shaihi@post.tau.ac.il>
//                 Eugene Lipovetsky <eug@post.tau.ac.il>
//                 Efi Fogel         <efif@post.tau.ac.il>
//                 Ron Wein          <wein@post.tau.ac.il>
#ifndef CGAL_PM_SEGMENT_DEBUG_TRAITS_2_H
#define CGAL_PM_SEGMENT_DEBUG_TRAITS_2_H

#include <CGAL/tags.h>

CGAL_BEGIN_NAMESPACE

template <class Kernel_>
class Pm_segment_traits_2 : public Kernel_190
{
  class My_x_curve_2{
  public:
    typedef typename Kernel_::Point_2     Point_2;
    typedef typename Kernel_::Segment_2   Segment;

    Segment seg;

    My_x_curve_2()
    {}

    My_x_curve_2(Segment s): seg(s)
    {}

    My_x_curve_2(Point_2 p1, Point_2 p2): seg(p1,p2)
    {}
  };

public:
  typedef Kernel_                         Kernel;

  // Categories:
  // #define HAS_LEFT_NOT
#if !defined(HAS_LEFT_NOT)
  typedef Tag_true                        Has_left_category;
#else
  typedef Tag_false                       Has_left_category;
#endif
    
  // #define HAS_REFLECT
#if !defined(HAS_REFLECT)
  typedef Tag_false                       Has_reflect_category;
#else
  typedef Tag_true                        Has_reflect_category;
#endif

  // Traits objects
  typedef typename Kernel::Point_2        Point_2;
  typedef My_x_curve_2                    X_monotone_curve_2;

  // Backward compatability    
  typedef Point_2                         Point;
  typedef X_monotone_curve_2              X_curve;

protected:
  // Functors:
  typedef typename Kernel::Is_vertical_2        Is_vertical_2;
  typedef typename Kernel::Construct_vertex_2   Construct_vertex_2;
  typedef typename Kernel::Less_x_2             Less_x_2;
  typedef typename Kernel::Equal_2              Equal_2;
    
public:
  // Creation
  Pm_segment_traits_2() {}

  // Operations
  // ----------
    
  /*! compare_x() compares the x-coordinates of two given points
   * \param p1 the first point
   * \param p2 the second point
   * \return LARGER if x(p1) > x(p2); SMALLER if x(p1) < x(p2); or else EQUAL
   */
  Comparison_result compare_x(const Point_2 & p1, const Point_2 & p2) const
  { return compare_x_2_object()(p1, p2); }

  /*! compare_xy() compares lexigoraphically the two points by x, then by y.
   * \param p1 the first point
   * \param p2 the second point
   * \return LARGER if x(p1) > x(p2), or if x(p1) = x(p2) and y(p1) > y(p2); 
   *         SMALLER if x(p1) < x(p2), or if x(p1) = x(p2) and y(p1) < y(p2);
   *         or else EQUAL
   */
  Comparison_result compare_xy(const Point_2 & p1, const Point_2 & p2) const
  { return compare_xy_2_object()(p1, p2); }

  /*! curve_is_vertical()
   * \param cv the curve
   * \return true iff the curve is vertical
   */
  bool curve_is_vertical(const X_monotone_curve_2 & cv) const 
  { return is_vertical_2_object()(cv.seg); }

  /*! point_in_x_range()
   * \param cv the curve
   * \param q the point
   * \return true if q is in the x range of cv
   *
   * \todo Intorduce Is_in_x_range_2() or perhaps Is_in_x_closed_range_2()
   * in kernel. Currently, this is implemented using existing traits (kernel)
   * functions (curve_source(), curve_target()) that return the source and
   * target points by value, which is not as efficient as possible.
   */
  bool point_in_x_range(const X_monotone_curve_2 & cv, const Point_2 & q) const
  {
    Construct_vertex_2 construct_vertex = construct_vertex_2_object();
    const Point_2 & source = construct_vertex(cv.seg, 0);
    const Point_2 & target = construct_vertex(cv.seg, 1);
    Less_x_2 less_x = less_x_2_object();
    return !((less_x(source, q) && less_x(target, q)) ||
             (less_x(q, source) && less_x(q, target)));
  }

  /*! curves_compare_y_at_x() compares the y-coordinate of two given curves at
   * the x-coordinate of a given point.
   * \param cv1 the first curve
   * \param cv2 the second curve
   * \param q the point
   * \return LARGER if cv1(x(q)) > cv2(x(q)); SMALLER if cv1(x(q)) < cv2(x(q));
   * or else EQUAL.
   * \pre The point q is in the x range of the two curves.
   */
  Comparison_result curves_compare_y_at_x(const X_monotone_curve_2 & cv1, 
                                          const X_monotone_curve_2 & cv2, 
                                          const Point_2 & q) const
  {
    CGAL_precondition(point_in_x_range(cv1.seg, q));
    CGAL_precondition(point_in_x_range(cv2.seg, q));

    return compare_y_at_x_2_object()(q, cv1.seg, cv2.seg);
  }

#if !defined(HAS_LEFT_NOT)
  /*! curves_compare_y_at_x_left() compares the y value of two curves in an
   * epsilon environment to the left of the x value of the input point
   * Preconditions: The point q is in the x range of the two curves, and both
   * of them must be also be defined to its left. The two curves must also
   * intersect at x(q).
   */
  Comparison_result curves_compare_y_at_x_left(const X_monotone_curve_2 & cv1,
                                               const X_monotone_curve_2 & cv2, 
                                               const Point_2 & q) const 
  {
    // The two curves must not be vertical.
    CGAL_precondition(! curve_is_vertical(cv1.seg));
    CGAL_precondition(! curve_is_vertical(cv2.seg));

    // The two curve must be defined at q and also to its left.
    CGAL_precondition_code(
        Construct_vertex_2 construct_vertex = construct_vertex_2_object();
	Less_x_2 less_x = less_x_2_object();
	const Point_2 & source1 = construct_vertex(cv1.seg, 0);
	const Point_2 & target1 = construct_vertex(cv1.seg, 1);
	const Point_2 & source2 = construct_vertex(cv2.seg, 0);
	const Point_2 & target2 = construct_vertex(cv2.seg, 1);
	);

    CGAL_precondition (less_x(source1, q) || less_x(target1, q));
    CGAL_precondition (!(less_x(source1, q) && less_x(target1, q)));
    
    CGAL_precondition (less_x(source2, q) || less_x(target2, q));
    CGAL_precondition (!(less_x(source2, q) && less_x(target2, q)));
    
    // Since the curves are continuous, if they are not equal at q, the same
    // result also applies to q's left.
    CGAL_precondition (compare_y_at_x_2_object()(q, cv1.seg, cv2.seg)==EQUAL);
    
    // <cv2> and <cv1> meet at a point with the same x-coordinate as q
    // compare their derivatives.
    return compare_slope_2_object()(cv2.seg, cv1.seg);
  }
#else
  /*! point_reflect_in_x_and_y() reflects the given point about the origin
   */
  Point_2 point_reflect_in_x_and_y(const Point_2 & pt) const
  {
    Point_2 org = construct_point_2_object()(ORIGIN);      
    typename Kernel::Vector_2 v = construct_vector_2_object()(pt, org);
    Point_2 reflected_pt(v);
    return reflected_pt;
  }

  /*! curve_reflect_in_x_and_y reflects the given curve about the origin
   */
  X_monotone_curve_2
  curve_reflect_in_x_and_y(const X_monotone_curve_2 & cv) const
  {
    X_monotone_curve_2 reflected_cv(point_reflect_in_x_and_y(cv.source()),
                                    point_reflect_in_x_and_y(cv.target()));
    return reflected_cv;
  }
#endif
    
  /*! curves_compare_y_at_x_right() compares the y value of two curves in an
   * epsilon environment to the right of the x value of the input point
   * Preconditions: The point q is in the x range of the two curves, and both
   * of them must be also be defined to its right. The two curves must also
   * intersect at x(q).
   */
  Comparison_result curves_compare_y_at_x_right(const X_monotone_curve_2 & cv1,
                                                const X_monotone_curve_2 & cv2,
                                                const Point_2 & q) const
  {
    // The two curves must not be vertical.
    CGAL_precondition(! curve_is_vertical(cv1));
    CGAL_precondition(! curve_is_vertical(cv2));

    // The two curve must be defined at q and also to its right.
    CGAL_precondition_code(
        Construct_vertex_2 construct_vertex = construct_vertex_2_object();
	Less_x_2 less_x = less_x_2_object();
	const Point_2 & source1 = construct_vertex(cv1.seg, 0);
	const Point_2 & target1 = construct_vertex(cv1.seg, 1);
	const Point_2 & source2 = construct_vertex(cv2.seg, 0);
	const Point_2 & target2 = construct_vertex(cv2.seg, 1);
	);

    CGAL_precondition (less_x(q, source1) || less_x(q, target1));
    CGAL_precondition (!(less_x(q, source1) && less_x(q, target1)));
    
    CGAL_precondition (less_x(q, source2) || less_x(q, target2));
    CGAL_precondition (!(less_x(q, source2) && less_x(q, target2)));
    
    // Since the curves are continuous, if they are not equal at q, the same
    // result also applies to q's left.
    CGAL_precondition (curves_compare_y_at_x(cv1, cv2, q) == EQUAL);     
    
    // <cv1> and <cv2> meet at a point with the same x-coordinate as q
    // compare their derivatives
    return compare_slope_2_object()(cv1.seg, cv2.seg);
  }
    
  /*! Return the location of the given point with respect to the input curve.
   * \param cv The curve.
   * \param p The point.
   * \pre p is in the x-range of cv.
   * \return LARGER if y(p) > cv(x(p));
   *         SMALLER if y(p) < cv(x(p));
   *         or else (if p is on the curve) EQUAL.
   */
  Comparison_result curve_compare_y_at_x (const Point_2 & p,
                                          const X_monotone_curve_2 & cv) const
  {
    CGAL_precondition(point_in_x_range(cv.seg, p));
    return compare_y_at_x_2_object()(p, cv.seg);

  }

  /*! Check if the two curves are the same (have the same graph).
   * \param cv1 The first curve.
   * \param cv2 The second curve.
   * \return (true) if the two curves are the same.
   */
  bool curve_equal(const X_monotone_curve_2 & cv1,
                   const X_monotone_curve_2 & cv2) const
  {
    Equal_2 equal = equal_2_object();
    const X_monotone_curve_2 & ocv1 = 
      construct_opposite_segment_2_object()(cv1.seg);
    return equal(cv1.seg, cv2.seg) || equal(ocv1.seg, cv2.seg);
  }

  /*! Check if the two points are the same.
   * \param p1 The first point.
   * \param p2 The second point.
   * \return (true) if p1 == p2.
   */
  bool point_equal(const Point_2 & p1, const Point_2 & p2) const
  { return equal_2_object()(p1, p2); }
  
  /*! Get the curve source.
   * \param cv The curve.
   * \return The source point.
   */
  Point_2 curve_source(const X_monotone_curve_2 & cv) const 
  { return construct_vertex_2_object()(cv.seg, 0); }

  /*! Get the curve target.
   * \param cv The curve.
   * \return The target point.
   */
  Point_2 curve_target(const X_monotone_curve_2 & cv) const 
  { return construct_vertex_2_object()(cv.seg, 1); }
};

CGAL_END_NAMESPACE

#endif
