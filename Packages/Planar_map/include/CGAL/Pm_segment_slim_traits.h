#ifndef CGAL_PM_SEGMENT_TRAITS_H
#define CGAL_PM_SEGMENT_TRAITS_H

#include <CGAL/Planar_map/Pm_segment_utilities_2.h>

CGAL_BEGIN_NAMESPACE

template <class Kernel_>
class Pm_segment_traits : public Kernel_
{
public:

  typedef Kernel_                      Kernel;

  // traits objects
  typedef typename Kernel::Point_2     Point_2;
  typedef Point_2                      Point; // for backward compatability

  typedef typename Kernel::Segment_2   X_curve;

  // Things I get from the kernel
  // ----------------------------

  //   Future interface:
  //  
  typedef typename Kernel::Is_vertical_2      Is_vertical_2;
  typedef typename Kernel::Equal_2            Equal_2;
  typedef typename Kernel::Has_on_2           Has_on_2;
  typedef typename Kernel::Compare_x_2        Compare_x_2;
  typedef typename Kernel::Compare_y_2        Compare_y_2;
  typedef typename Kernel::Construct_vertex_2 Construct_vertex_2;
  typedef typename Kernel::Construct_opposite_direction_2
                                              Construct_opposite_direction_2;
  typedef CGAL::Construct_direction_at_endpoint_2<Kernel, X_curve>
                                              Construct_direction_2;
  typedef CGAL::Compare_y_at_x_for_segments_2<Kernel, X_curve>
                                              Compare_y_at_x_2;
  typedef CGAL::Counterclockwise_in_between_for_segments_2<Kernel, X_curve>
                                              Counterclockwise_in_between_2;
  

  inline Construct_direction_2 construct_direction_2_object() const
  { return Construct_direction_2(); }

  inline Compare_y_at_x_2 compare_y_at_x_2_object() const
  { return Compare_y_at_x_2(); }

  inline Counterclockwise_in_between_2 counterclockwise_in_between_2_object() const
  { return Counterclockwise_in_between_2(); }

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

public:
  // Creation
  Pm_segment_traits() {}

  // SHAI: What about kernels with state?
  
  // Compare the y value of two curves in an epsilon environment to
  // the left of the x value of the input point
  Comparison_result 
  curve_compare_at_x_left(const X_curve &cv1, const X_curve &cv2, 
                          const Point_2 &q) const 
  {
    // If one of the curves is vertical then return EQUAL.
    if ( is_vertical_2_object()(cv1) || 
	 (is_vertical_2_object()(cv2)) ) 
      return EQUAL;
    // If one of the curves is not defined at q then return EQUAL.
    if ( ! is_left(leftmost(cv1.source(), cv1.target()), q) ) return EQUAL;
    if ( ! is_left(leftmost(cv2.source(), cv2.target()), q) ) return EQUAL;
    
    Comparison_result r = compare_y_at_x_2_object()(q, cv1, cv2);
        
    if ( r != EQUAL )
      return r;     // since the curve is continous 
    
    // <cv2> and <cv1> meet at a point with the same x-coordinate as q
    // compare their derivatives
    return compare_slope_2_object()(cv2, cv1);
  }
  
  // Compare the y value of two curves in an epsilon environment to
  // the right of the x value of the input point
  Comparison_result 
  curve_compare_at_x_right(const X_curve &cv1, const X_curve &cv2, 
			   const Point_2 & q) const 
  {
    // If one of the curves is vertical then return EQUAL.
    if ( is_vertical_2_object()(cv1) || 
	 (is_vertical_2_object()(cv2))  ) 
      return EQUAL;
    // If one of the curves is not defined at q then return EQUAL.
    if ( ! is_right(rightmost(cv1.source(), cv1.target()), q) ) return EQUAL;
    if ( ! is_right(rightmost(cv2.source(), cv2.target()), q) ) return EQUAL;
    
    Comparison_result r = compare_y_at_x_2_object()(q, cv1, cv2);
    
    if ( r != EQUAL)
      return r;     // since the curve is continous 
    
    // <cv1> and <cv2> meet at a point with the same x-coordinate as q
    // compare their derivatives
    return compare_slope_2_object()(cv1, cv2);
  }
 
};

CGAL_END_NAMESPACE

#endif // CGAL_PM_SEGMENT_EXACT_TRAITS_H
