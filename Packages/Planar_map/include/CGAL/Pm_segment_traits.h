// ======================================================================
//
// Copyright (c) 2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.4-I-40 $
// release_date  : $CGAL_Date: 2001/12/28 $
//
// file          : include/CGAL/Pm_segment_traits.h
// package       : Planar_map (5.80)
// maintainer    : Eyal Flato <flato@math.tau.ac.il>
// author(s)     : Oren Nechushtan <theoren@math.tau.ac.il>
//                 Iddo Hanniel <hanniel@math.tau.ac.il>
//                 Shai Hirsch <shaihi@post.tau.ac.il>
//
// coordinator   : Tel-Aviv University (Dan Halperin halperin<@math.tau.ac.il>)
//
// ======================================================================
#ifndef CGAL_PM_SEGMENT_TRAITS_H
#define CGAL_PM_SEGMENT_TRAITS_H

// Status on Dec. 4th, 2001
// Class was converted to use as much of the kernel as currently possible

CGAL_BEGIN_NAMESPACE


template <class Kernel_>
class Pm_segment_traits
{
public:

  typedef Kernel_                    Kernel;

  // traits objects
  typedef typename Kernel::Point_2   Point_2;
  typedef Point_2                    Point; // for backward compatability

  typedef typename Kernel::Segment_2 X_curve;

  // Things I get from the kernel
  // ----------------------------
  //
  //   Future interface:
  //  
  //   typedef typename Kernel::Is_vertical_2      Is_vertical_2;
  //   typedef typename Kernel::Compare_y_at_x_2   Compare_y_at_x_2;
  //   typedef typename Kernel::Counterclockwise_in_between_2
  //                                               Counterclockwise_in_between_2;
  //   typedef typename Kernel::Equal_2            Equal_2;
  //   typedef typename Kernel::Has_on_2           Has_on_2;
  //   typedef typename Kernel::Compare_x_2        Compare_x_2;
  //   typedef typename Kernel::Compare_y_2        Compare_y_2;
  //   typedef typename Kernel::Construct_vertex_2 Construct_vertex_2;


  //   Implementation
  //   
  //   typedef typename Kernel::Less_x_2                              Less_x_2;
  //   typedef typename Kernel::Construct_opposite_segment_2 
  //                                              Construct_opposite_segment_2;
  //   typedef typename Kernel::Construct_line_2              Construct_line_2;
  //   typedef typename Kernel::Construct_direction_2    Construct_direction_2;
  //   typedef typename Kernel::Construct_vertical_projected_point_2 
  //                                       onstruct_vertical_projected_point_2;

  // Currently, I leave this in the traits
  // Maybe we can change the usage inside Planar_map_2
  typedef enum
  {
    UNDER_CURVE        = -1,
    CURVE_NOT_IN_RANGE =  0,
    ABOVE_CURVE        =  1,
    ON_CURVE           =  2

  } Curve_point_status;	

private:
  Kernel m_kernel;

public:
  // Creation
  Pm_segment_traits() {}

  Pm_segment_traits(const Kernel& kernel) : m_kernel(kernel) {}
  
  // Access to curve source
  Point_2 curve_source(const X_curve & cv) const 
  { 
    return m_kernel.construct_vertex_2_object()(cv, 0);
  }
  
  // Access to curve target
  Point_2 curve_target(const X_curve & cv) const 
  {
    return m_kernel.construct_vertex_2_object()(cv, 1);
  }
  
  // Answers true iff the curve is vertical.
  bool curve_is_vertical(const X_curve & cv) const 
  {
    return m_kernel.is_vertical_2_object()(cv);
  }

  // Returns the curve-point status of the input objects
  Curve_point_status 
  curve_get_point_status(const X_curve &cv, const Point_2 & p) const
  {
    if ( ! curve_is_in_x_range(cv, p))
      return CURVE_NOT_IN_RANGE;
    if ( ! curve_is_vertical(cv))
      {
	// Calculate vertical projection on curve
	const Point_2 & proj = 
	  construct_vertical_projected_point_2_object(cv, p);
	int res = m_kernel.compare_y_2_object()(p, proj);
	if (res == SMALLER) return UNDER_CURVE;
	if (res == LARGER)	return ABOVE_CURVE;

	return ON_CURVE;
      }
    else
      {
	if (is_lower(p,lowest(curve_source(cv),curve_target(cv))))
	  return UNDER_CURVE;
	if (is_higher(p,highest(curve_source(cv),curve_target(cv))))
	  return ABOVE_CURVE;

	return ON_CURVE;
      }
  }

  // Compares the y value of two curves at an x value of the input point
  Comparison_result 
  curve_compare_at_x(const X_curve &cv1, const X_curve &cv2, const Point_2 &q) 
    const 
  {
    typename Kernel::Line_2 l1 = m_kernel.construct_line_2_object()(cv1);
    typename Kernel::Line_2 l2 = m_kernel.construct_line_2_object()(cv2);
    return m_kernel.compare_y_at_x_2_object()(q, l1, l2);
  }

  // Compare the y value of two curves in an epsilon environment to
  // the left of the x value of the input point
  Comparison_result 
  curve_compare_at_x_left(const X_curve &cv1, const X_curve &cv2, 
                          const Point_2 &q) const 
  {
    // If one of the curves is vertical then return EQUAL.
    if ( curve_is_vertical(cv1) || (curve_is_vertical(cv2)) ) return EQUAL;
    // If one of the curves is not defined at q then return EQUAL.
    if ( ! is_left(leftmost(cv1.source(), cv1.target()), q) ) return EQUAL;
    if ( ! is_left(leftmost(cv2.source(), cv2.target()), q) ) return EQUAL;
    
    Comparison_result r = curve_compare_at_x(cv1, cv2, q);
        
    if ( r != EQUAL )
      return r;     // since the curve is continous 
    
    // <cv2> and <cv1> meet at a point with the same x-coordinate as q
    // compare their derivatives
    return m_kernel.compare_slope_2_object()(cv2, cv1);
    // return compare_value(curve_derivative(cv2), curve_derivative(cv1));
  }
  
  // Compare the y value of two curves in an epsilon environment to
  // the right of the x value of the input point
  Comparison_result 
  curve_compare_at_x_right(const X_curve &cv1, const X_curve &cv2, 
			   const Point_2 & q) const 
  {
    // If one of the curves is vertical then return EQUAL.
    if ( curve_is_vertical(cv1) || (curve_is_vertical(cv2))  ) return EQUAL;
    // If one of the curves is not defined at q then return EQUAL.
    if ( ! is_right(rightmost(cv1.source(), cv1.target()), q) ) return EQUAL;
    if ( ! is_right(rightmost(cv2.source(), cv2.target()), q) ) return EQUAL;
    
    Comparison_result r = curve_compare_at_x(cv1, cv2, q);
    
    if ( r != EQUAL)
      return r;     // since the curve is continous 
    
    // <cv1> and <cv2> meet at a point with the same x-coordinate as q
    // compare their derivatives
    return m_kernel.compare_slope_2_object()(cv1, cv2);
    //return compare_value(curve_derivative(cv1), curve_derivative(cv2));
  }
  
  // done
  // (Check what happens if cv == first, if first == second 
  // and if both.)
  bool curve_is_between_cw(const X_curve &cv, 
                           const X_curve &first, 
                           const X_curve &second, 
                           const Point_2 &point) const
  {
    typedef typename Kernel::Direction_2 Direction_2;
    
    X_curve my_cv = cv, my_first = first, my_second = second;
    if ( curve_source(my_cv)    != point ) my_cv     = curve_flip(cv);
    if ( curve_source(my_first) != point ) my_first  = curve_flip(first);
    if ( curve_source(my_second)!= point ) my_second = curve_flip(second);

    Direction_2 d  = m_kernel.construct_direction_2_object()(my_cv);
    Direction_2 d1 = m_kernel.construct_direction_2_object()(my_first);
    Direction_2 d2 = m_kernel.construct_direction_2_object()(my_second);

    return m_kernel.counterclockwise_in_between_2_object()(d, d1, d2);
  }

  // Compares the x value of two points
  Comparison_result compare_x(const Point_2 &p1, const Point_2 &p2) const
  { return m_kernel.compare_x_2_object()(p1, p2); }

  // Compares the y value of two points
  Comparison_result compare_y(const Point_2 &p1, const Point_2 &p2) const
  { return m_kernel.compare_y_2_object()(p1, p2); }

  bool curve_is_same(const X_curve & cv1,const X_curve & cv2) const
  { return m_kernel.equal_2_object()(cv1, cv2); }

  // Intorduce Is_in_x_range_2 / Is_in_x_closed_range_2 ?
  // This can be implemented on the traits_wrap level by using other
  // simpler predicated from the Kernel / traits.
  // Used in the Bounding Box, but there it probably get it from
  // the Pm's Traits_wrap
  bool curve_is_in_x_range(const X_curve & cv, const Point_2 & q) const
  { 
    return !( is_right(q, rightmost(cv.source(), cv.target())) ||
              is_left(q, leftmost(cv.source(), cv.target()))	 );
  }

private:
  // constructs the opposite segment (with the source and target
  // exchanged)
  // Used internally and in the Arrangement, so shouldn't be part of
  // this interface
  X_curve curve_flip(const X_curve &cv) const
  {
    return m_kernel.construct_opposite_segment_2_object()(cv);
  }
  // These stuff need to be cached
  bool is_left(const Point_2 &p1, const Point_2 &p2) const 
  { return m_kernel.less_x_2_object()(p1, p2); }
  bool is_right(const Point_2 &p1, const Point_2 &p2) const 
  { return m_kernel.less_x_2_object()(p2, p1); }
  bool is_same_x(const Point_2 &p1, const Point_2 &p2) const 
  { return m_kernel.equal_x_object()(p1, p2); }
  bool is_lower(const Point_2 &p1, const Point_2 &p2) const 
  { return m_kernel.less_y_2_object()(p1, p2); }
  bool is_higher(const Point_2 &p1, const Point_2 &p2) const 
  { return m_kernel.less_y_2_object()(p2, p1); }
  bool is_same_y(const Point_2 &p1, const Point_2 &p2) const 
  { return m_kernel.equal_y_object()(p1, p2); }
  bool is_same(const Point_2 &p1, const Point_2 &p2) const
  {
    return (compare_x(p1, p2) == EQUAL) &&
      (compare_y(p1, p2) == EQUAL);
  }
  const Point_2& leftmost(const Point_2 &p1, const Point_2 &p2) const
  { return (is_left(p1, p2) ? p1 : p2); }

  const Point_2& rightmost(const Point_2 &p1, const Point_2 &p2) const
  { return (is_right(p1, p2) ? p1 : p2); }
  
  const Point_2& lowest(const Point_2 &p1, const Point_2 &p2) const
  { return (is_lower(p1, p2) ? p1 : p2); }
  
  const Point_2& highest(const Point_2 &p1, const Point_2 &p2) const
  { return (is_higher(p1, p2) ? p1 : p2); }
  
  // Comment this one ! ##############
  
  Point_2 construct_vertical_projected_point_2_object(const X_curve &cv, const Point_2 & q) const
    {
      if ( ! curve_is_in_x_range(cv, q) )
	return cv.source();
      
      if (curve_is_vertical(cv))
	return cv.source();
      
      const Point_2 & a = cv.source();
      const Point_2 & b = cv.target();
      return Point_2 ((b.hx() * a.hw() - a.hx() * b.hw()) * q.hx() * a.hw(),
		      (b.hx() * a.hw() - a.hx() * b.hw()) * q.hw() * a.hy() + 
		      (b.hy() * a.hw() - a.hy() * b.hw()) * 
		      (q.hx() * a.hw() - a.hx() * q.hw()),  
		      (b.hx() * a.hw() - a.hx() * b.hw()) * q.hw() * a.hw());
    }

};

CGAL_END_NAMESPACE

#endif // CGAL_PM_SEGMENT_EXACT_TRAITS_H
