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
// author(s)     : Iddo Hanniel      <hanniel@math.tau.ac.il>
//                 Eyal Flato        <flato@post.tau.ac.il>
//                 Oren Nechushtan   <theoren@math.tau.ac.il>
//                 Eti Ezra          <estere@post.tau.ac.il>
//                 Shai Hirsch       <shaihi@post.tau.ac.il>
//                 Eugene Lipovetsky <eug@post.tau.ac.il>
//                 Efi Fogel         <efif@post.tau.ac.il>
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ======================================================================
#ifndef CGAL_PM_SEGMENT_TRAITS_2_H
#define CGAL_PM_SEGMENT_TRAITS_2_H

#include <CGAL/Planar_map_2/Pm_segment_utilities_2.h>
#include <CGAL/tags.h>

CGAL_BEGIN_NAMESPACE

template <class Kernel_>
class Pm_segment_traits_2 : public Kernel_
{
public:
  typedef Kernel_                         Kernel;

  // Categories:
#define HAS_LEFT_NOT
#if !defined(HAS_LEFT_NOT)
  typedef Tag_true                        Has_left_category;
#else
  typedef Tag_false                       Has_left_category;
#endif
    
  // Traits objects
  typedef typename Kernel::Point_2        Point_2;
  typedef typename Kernel::Segment_2      X_curve_2;

  // Backward compatability    
  typedef Point_2                         Point;
  typedef X_curve_2                       X_curve;

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
    
  typedef CGAL::Counterclockwise_in_between_for_segments_2<Kernel, X_curve_2>
                                                Counterclockwise_in_between_2;

protected:
  inline Counterclockwise_in_between_2 counterclockwise_in_between_2_object()
    const
  { return Counterclockwise_in_between_2(); }

public:
  // Creation
  Pm_segment_traits_2() {}

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
  { return compare_x_2_object()(p1, p2); }

  /*! compare_y() compares the y-coordinates of two given points
   * \param p1 the first point
   * \param p2 the second point
   * \return LARGER if y(p1) > y(p2), SMALLER if y(p1) < y(p2), or else EQUAL
   *
   * \todo replace indirect use compare_y() with compare_y_2()
   */
  Comparison_result compare_y(const Point_2 & p1, const Point_2 & p2) const
  { return compare_y_2_object()(p1, p2); }

  /*! curve_is_vertical()
   * \param cv the curve
   * \return true iff the curve is vertical
   *
   * \todo replace indirect use curve_is_vertical() with is_vertical_2()
   */
  bool curve_is_vertical(const X_curve_2 & cv) const 
  { return is_vertical_2_object()(cv); }

  /*! curve_is_in_x_range()
   * \param cv the curve
   * \param q the point
   * \return true if q is in the x range of cv
   *
   * \todo Intorduce Is_in_x_range_2() or perhaps Is_in_x_closed_range_2()
   * in kernel. Currently, this is implemented using existing traits (kernel)
   * functions (curve_source(), curve_target()) that return the source and
   * target points by value, which is not as efficient as possible.
   */
  bool curve_is_in_x_range(const X_curve_2 & cv, const Point_2 & q) const
  {
    Construct_vertex_2 construct_vertex = construct_vertex_2_object();
    const Point_2 & source = construct_vertex(cv, 0);
    const Point_2 & target = construct_vertex(cv, 1);
    Less_x_2 less_x = less_x_2_object();
    return !((less_x(source, q) && less_x(target, q)) ||
             (less_x(q, source) && less_x(q, target)));
  }

  /*! curve_compare_at_x() compares the y-coordinate of two given curves at
   * the x-coordinate of a given point
   * \param cv1 the first curve
   * \param cv2 the second curve
   * \param q the point
   * \return EQUAL if at least one of cv1 and cv2 is not defined at q's
   * x-coordinate x(q). Otherwise, LARGER if cv1(x(q)) > cv2(x(q)), SMALLER if
   * cv1(x(q)) < cv2(x(q), or else EQUAL.
   * \todo replace indirect use curve_compare_at_x() with compare_y_at_x_2()
   */
  Comparison_result curve_compare_at_x(const X_curve_2 & cv1, 
				       const X_curve_2 & cv2, 
				       const Point_2 & q) const
  {
    if (!curve_is_in_x_range(cv1, q) || !curve_is_in_x_range(cv2, q))
      return EQUAL;
    return compare_y_at_x_2_object()(q, cv1, cv2);
  }

#if !defined(HAS_LEFT_NOT)
  /*! curve_compare_at_x_left() compares the y value of two curves in an
   * epsilon environment to the left of the x value of the input point
   */
  Comparison_result curve_compare_at_x_left(const X_curve_2 & cv1,
                                            const X_curve_2 & cv2, 
                                            const Point_2 & q) const 
  {
    // If one of the curves is vertical then return EQUAL.
    Is_vertical_2 is_vertical = is_vertical_2_object();
    if (is_vertical(cv1) || (is_vertical(cv2))) return EQUAL;

    // If one of the curves is not defined at q then return EQUAL.
    Construct_vertex_2 construct_vertex = construct_vertex_2_object();
    Less_x_2 less_x = less_x_2_object();
    const Point_2 & source1 = construct_vertex(cv1, 0);
    const Point_2 & target1 = construct_vertex(cv1, 1);
    if (!(less_x(source1, q) || less_x(target1, q))) return EQUAL;
    
    const Point_2 & source2 = construct_vertex(cv2, 0);
    const Point_2 & target2 = construct_vertex(cv2, 1);
    if (!(less_x(source2, q) || less_x(target2, q))) return EQUAL;

    if (less_x(source1, q) && less_x(target1, q)) return EQUAL;
    if (less_x(source2, q) && less_x(target2, q)) return EQUAL;
    
    // since the curve is continous 
    Comparison_result r = compare_y_at_x_2_object()(q, cv1, cv2);
    if (r != EQUAL) return r;     
    
    // <cv2> and <cv1> meet at a point with the same x-coordinate as q
    // compare their derivatives
    return compare_slope_2_object()(cv2, cv1);
  }
#else
  /*! point_reflect_in_x_and_y() reflects the given point about the origin
   */
  Point_2 point_reflect_in_x_and_y(const Point_2 & pt) const
  {
    Point_2 org = construct_point_2_object()(ORIGIN);      
    typename Kernel::Vector_2 v = construct_vector_2_object()(org, pt);
    Point_2 reflected_pt(v);
    return reflected_pt;
  }

  /*! curve_reflect_in_x_and_y reflects the given curve about the origin
   */
  X_curve_2 curve_reflect_in_x_and_y(const X_curve_2 & cv) const
  {
    X_curve_2 reflected_cv(point_reflect_in_x_and_y ( cv.source()),
                           point_reflect_in_x_and_y ( cv.target()));
    return reflected_cv;
  }
#endif
    
  /*! curve_compare_at_x_right() compares the y value of two curves in an
   * epsilon environment to the right of the x value of the input point
   */
  Comparison_result curve_compare_at_x_right(const X_curve_2 & cv1,
                                             const X_curve_2 & cv2, 
                                             const Point_2 & q) const 
  {
    // If one of the curves is vertical then return EQUAL.
    Is_vertical_2 is_vertical = is_vertical_2_object();
    if (is_vertical(cv1) || (is_vertical(cv2))) return EQUAL;

    // If one of the curves is not defined at q then return EQUAL.
    Construct_vertex_2 construct_vertex = construct_vertex_2_object();
    Less_x_2 less_x = less_x_2_object();
    const Point_2 & source1 = construct_vertex(cv1, 0);
    const Point_2 & target1 = construct_vertex(cv1, 1);
    if (!(less_x(q, source1) || less_x(q, target1))) return EQUAL;

    const Point_2 & source2 = construct_vertex(cv2, 0);
    const Point_2 & target2 = construct_vertex(cv2, 1);
    if (!(less_x(q, source2) || less_x(q, target2))) return EQUAL;

    if (less_x(q, source1) && less_x(q, target1)) return EQUAL;
    if (less_x(q, source2) && less_x(q, target2)) return EQUAL;
    
    // since the curve is continous (?)
    Comparison_result r = curve_compare_at_x(cv1, cv2, q);
    if (r != EQUAL) return r;     
    
    // <cv1> and <cv2> meet at a point with the same x-coordinate as q
    // compare their derivatives
    return compare_slope_2_object()(cv1, cv2);
  }
    
  /*! Return the curve-point status of the input objects
   * \todo Remove curve_get_point_status() from wrapper. Verify that
   * curve_is_in_x_range() and compare_y_at_x_2() are required in the
   * traits, and use directly.
   */
  Curve_point_status 
  curve_get_point_status(const X_curve_2 & cv, const Point_2 & p) const
  {
    if (!curve_is_in_x_range(cv, p)) return CURVE_NOT_IN_RANGE;
    Comparison_result res = compare_y_at_x_2_object()(p, cv);
    return ((res == LARGER) ? ABOVE_CURVE :
            ((res == SMALLER) ? UNDER_CURVE : ON_CURVE));
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
    // Notice the change in order of first and second
    return counterclockwise_in_between_2_object()(point, cv, second, first);
  }

  /*! \todo replace indirect use curve_is_same() with equal_2()
   */
  bool curve_is_same(const X_curve_2 & cv1,const X_curve_2 & cv2) const
  { return equal_2_object()(cv1, cv2); }

  /*! \todo replace indirect use point_is_same() with equal_2()
   */
  bool point_is_same(const Point_2 & p1, const Point_2 & p2) const
  { return equal_2_object()(p1, p2); }
  
  /*! \todo replace indirect use curve_source() with construct_vertex_2()
   */
  Point_2 curve_source(const X_curve_2 & cv) const 
  { return construct_vertex_2_object()(cv, 0); }

  /*! \todo replace indirect use curve_source() with construct_vertex_2()
   */
  Point_2 curve_target(const X_curve_2 & cv) const 
  { return construct_vertex_2_object()(cv, 1); }
};

CGAL_END_NAMESPACE

#endif
