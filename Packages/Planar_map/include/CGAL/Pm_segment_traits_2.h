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
#ifndef CGAL_PM_SEGMENT_TRAITS_2_H
#define CGAL_PM_SEGMENT_TRAITS_2_H

#include <CGAL/tags.h>

CGAL_BEGIN_NAMESPACE

template <class Kernel_>
class Pm_segment_traits_2 : public Kernel_
{
public:
  typedef Kernel_                         Kernel;

  // Categories:
  //#define HAS_LEFT_NOT
#if !defined(HAS_LEFT_NOT)
  typedef Tag_true                        Has_left_category;
#else
  typedef Tag_false                       Has_left_category;
#endif

  //#define HAS_REFLECT
#if !defined(HAS_REFLECT)
  typedef Tag_false                       Has_reflect_category;
#else
  typedef Tag_true                        Has_reflect_category;
#endif
    
  // Traits objects
  typedef typename Kernel::Point_2        Point_2;
  typedef typename Kernel::Segment_2      X_monotone_curve_2;

  // Backward compatability    
  typedef Point_2                         Point;
  typedef X_monotone_curve_2              X_curve;

protected:
  // Functors:
  typedef typename Kernel::Compare_x_2          Compare_x_2;
  typedef typename Kernel::Compare_xy_2         Compare_xy_2;
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
  { return is_vertical_2_object()(cv); }

  /*! point_in_x_range()
   * \param cv the curve
   * \param q the point
   * \return true if q is in the x range of cv
   */
  bool point_in_x_range(const X_monotone_curve_2 & cv, const Point_2 & q) const
  {
#if 1
    Compare_x_2 compare_x = compare_x_2_object();
    Construct_vertex_2 construct_vertex = construct_vertex_2_object();
    Comparison_result res1 = compare_x(construct_vertex(cv, 0), q);
    Comparison_result res2 = compare_x(construct_vertex(cv, 1), q);

    // We check if x(p) equals the x value of one of the end-points.
    // If not, we check whether one end-point is to p's left and the other is
    // to its right.
    return ((res1 == EQUAL) || (res2 == EQUAL) || (res1 != res2));
#else
    // \todo use this code instead, after the calls it uses are supported
    // in the LEDA kernel.
    Compare_x_2 compare_x = compare_x_2_object();
    Comparison_result res1 = compare_x(cv, 0, q);
    Comparison_result res2 = compare_x(cv, 1, q);

    // We check if x(p) equals the x value of one of the end-points.
    // If not, we check whether one end-point is to p's left and the other is
    // to its right.
    return ((res1 == EQUAL) || (res2 == EQUAL) || (res1 != res2));
#endif
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
    CGAL_precondition(point_in_x_range(cv1, q));
    CGAL_precondition(point_in_x_range(cv2, q));

    return compare_y_at_x_2_object()(q, cv1, cv2);
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
    // The two curve must be defined at q and also to its left.
    CGAL_precondition_code(
        Construct_vertex_2 construct_vertex = construct_vertex_2_object();
	Compare_xy_2 cmp_xy = compare_xy_2_object();
	const Point_2 & source1 = construct_vertex(cv1, 0);
	const Point_2 & target1 = construct_vertex(cv1, 1);
	const Point_2 & source2 = construct_vertex(cv2, 0);
	const Point_2 & target2 = construct_vertex(cv2, 1);
        const Is_vertical_2 is_vertical = is_vertical_2_object();
	);

    CGAL_precondition((cmp_xy(source1, q) == SMALLER) ||
                      (cmp_xy(target1, q) == SMALLER));
    CGAL_precondition((cmp_xy(source1, q) != SMALLER) ||
                      (cmp_xy(target1, q) != SMALLER));
    
    CGAL_precondition((cmp_xy(source2, q) == SMALLER) ||
                      (cmp_xy(target2, q) == SMALLER));
    CGAL_precondition((cmp_xy(source2, q) != SMALLER) ||
                      (cmp_xy(target2, q) != SMALLER));
    
    // Since the curves are continuous, if they are not equal at q, the same
    // result also applies to q's left.
    CGAL_precondition((is_vertical(cv1) && has_on_2_object()(cv1, q)) ||
                      (is_vertical(cv2) && has_on_2_object()(cv2, q)) ||
                      (compare_y_at_x_2_object()(q, cv1, cv2) == EQUAL));
    
    // <cv2> and <cv1> meet at a point with the same x-coordinate as q
    // compare their derivatives.
    return compare_slope_2_object()(cv2, cv1);
  }
#else
  /*! point_reflect_in_x_and_y() reflects the given point about the origin
   */
  Point_2 point_reflect_in_x_and_y(const Point_2 & pt) const
  {
    Point_2 org = construct_point_2_object()(ORIGIN);      
    typename Kernel::Vector_2 v = construct_vector_2_object()(pt, org);
    Point_2 reflected_pt = org + v;
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
    // The two curve must be defined at q and also to its right.
    CGAL_precondition_code(
        Construct_vertex_2 construct_vertex = construct_vertex_2_object();
	Compare_xy_2 cmp_xy = compare_xy_2_object();
	const Point_2 & source1 = construct_vertex(cv1, 0);
	const Point_2 & target1 = construct_vertex(cv1, 1);
	const Point_2 & source2 = construct_vertex(cv2, 0);
	const Point_2 & target2 = construct_vertex(cv2, 1);
	);

    CGAL_precondition((cmp_xy(q, source1) == SMALLER)||
                      (cmp_xy(q, target1) == SMALLER));
    CGAL_precondition((cmp_xy(q, source1) != SMALLER) |
                      (cmp_xy(q, target1) != SMALLER));
    
    CGAL_precondition((cmp_xy(q, source2) == SMALLER) ||
                      (cmp_xy(q, target2) == SMALLER));
    CGAL_precondition((cmp_xy(q, source2) != SMALLER) ||
                      (cmp_xy(q, target2) != SMALLER));
    
    // Since the curves are continuous, if they are not equal at q, the same
    // result also applies to q's left.
    CGAL_precondition((is_vertical_2_object()(cv1) &&
                       has_on_2_object()(cv1, q)) ||
                      (is_vertical_2_object()(cv2) &&
                       has_on_2_object()(cv2, q)) ||
                      (compare_y_at_x_2_object()(q, cv1, cv2) == EQUAL));
    
    // <cv1> and <cv2> meet at a point with the same x-coordinate as q
    // compare their derivatives
    return compare_slope_2_object()(cv1, cv2);
  }
    
  /*! Return the location of the given point with respect to the input curve.
   * \param cv The curve.
   * \param p The point.
   * \pre p is in the x-range of cv.
   * \return SMALLER if y(p) < cv(x(p)), that is the point is below the curve;
   *         LARGER if y(p) > cv(x(p)), that is the point is above the curve;
   *         or else (if p is on the curve) EQUAL.
   */
  Comparison_result curve_compare_y_at_x (const Point_2 & p,
                                          const X_monotone_curve_2 & cv) const
  {
    CGAL_precondition(point_in_x_range(cv, p));
    return compare_y_at_x_2_object()(p, cv);
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
      construct_opposite_segment_2_object()(cv1);
    return equal(cv1, cv2) || equal(ocv1, cv2);
  }

  /*! Check if the two points are the same.
   * \param p1 The first point.
   * \param p2 The second point.
   * \return (true) if p1 == p2.
   */
  bool point_equal(const Point_2 & p1, const Point_2 & p2) const
  { return equal_2_object()(p1, p2); }
  
  /*! Obtain the curve source.
   * We return the point by value (and by reference), because the implementation
   * of the Construct_vertex_2 function object in undelying kernel may return
   * a temporary variable.
   * \param cv The curve.
   * \return The source point.
   */
  const Point_2 curve_source(const X_monotone_curve_2 & cv) const 
  { return construct_vertex_2_object()(cv, 0); }

  /*! Obtain the curve target.
   * We return the point by value (and by reference), because the implementation
   * of the Construct_vertex_2 function object in undelying kernel may return
   * a temporary variable.
   * \param cv The curve.
   * \return The target point.
   */
  const Point_2 curve_target(const X_monotone_curve_2 & cv) const 
  { return construct_vertex_2_object()(cv, 1); }
};

CGAL_END_NAMESPACE

#endif
