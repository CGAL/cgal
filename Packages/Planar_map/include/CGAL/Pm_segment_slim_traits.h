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
// file          : include/CGAL/Pm_segment_slim_traits.h
// package       : Planar_map (5.80)
// maintainer    : Eyal Flato <flato@math.tau.ac.il>
// author(s)     : Oren Nechushtan <theoren@math.tau.ac.il>
//                 Iddo Hanniel <hanniel@math.tau.ac.il>
//                 Shai Hirsch <shaihi@post.tau.ac.il>
//                 Efi Fogel <efifogel@post.tau.ac.il>
//
// coordinator   : Tel-Aviv University (Dan Halperin halperin<@math.tau.ac.il>)
//
// ======================================================================
#ifndef CGAL_PM_SEGMENT_TRAITS_H
#define CGAL_PM_SEGMENT_TRAITS_H

// Status on Dec. 4th, 2001
// Class was converted to use as much of the kernel as currently possible

CGAL_BEGIN_NAMESPACE

template <class Kernel_, class X_curve_>
struct Construct_direction_at_endpoint_2 : 
public Kernel_::Construct_direction_2
{
  typedef Kernel_                                Kernel;
  typedef X_curve_                               X_curve;

  typedef typename Kernel::Construct_direction_2 Base;
  typedef typename Kernel::Direction_2           Direction_2;
  typedef typename Kernel::Point_2               Point_2;

  Direction_2 operator()(const X_curve & cv, const Point_2 &) const
    { 
      return Base::operator()(cv);
    }
};

// Compares the y value of two curves at an x value of the input point
template <class Kernel_, class X_curve_>
struct Compare_y_at_x_for_segments_2 :
public Kernel_::Compare_y_at_x_2
{
  typedef Kernel_                           Kernel;
  typedef X_curve_                          X_curve;

  typedef typename Kernel::Compare_y_at_x_2 Base;
  typedef typename Kernel::Point_2          Point_2;
  typedef typename Kernel::Construct_line_2 Construct_line_2;

  Comparison_result 
  operator()( const Point_2 & q, const X_curve &cv1, const X_curve &cv2) const 
    {
      typename Kernel::Line_2 l1 = construct_line(cv1);
      typename Kernel::Line_2 l2 = construct_line(cv2);
      return Base::operator()(q, l1, l2);
    }

  Comparison_result 
  operator()( const Point_2 & q, const X_curve &cv) const 
    {
      /*
      if ( is_vertical_2_object()(cv) )
	{
	  const Point & src = construct_vertex_2_object()(cv, 0);
	  const Point & trg = construct_vertex_2_object()(cv, 1);
	  
	}
      */
      typename Kernel::Line_2 l = construct_line(cv);
      return Base::operator()(q, l);
    }
  
  Construct_line_2 construct_line;

};

/*
template <class Kernel>
inline
Pm_segment_traits<Kernel>::
Construct_direction_2 
Pm_segment_traits<Kernel>::
construct_direction_2_object() const
{
  return Construct_direction_2();
}
*/

template <class Kernel_>
class Pm_segment_traits : public Kernel_
{
public:

  typedef Kernel_                    Kernel;

  // traits objects
  typedef typename Kernel::Point_2     Point_2;
  typedef Point_2                      Point; // for backward compatability

  typedef typename Kernel::Direction_2 Direction_2;

  typedef typename Kernel::Segment_2   X_curve;

  // Things I get from the kernel
  // ----------------------------
  //
  //   Future interface:
  //  
  typedef typename Kernel::Is_vertical_2      Is_vertical_2;
  typedef typename Kernel::Counterclockwise_in_between_2
                                              Counterclockwise_in_between_2;
  typedef typename Kernel::Equal_2            Equal_2;
  typedef typename Kernel::Has_on_2           Has_on_2;
  typedef typename Kernel::Compare_x_2        Compare_x_2;
  typedef typename Kernel::Compare_y_2        Compare_y_2;
  typedef typename Kernel::Construct_vertex_2 Construct_vertex_2;
  
  typedef typename Kernel::Construct_opposite_direction_2
                                              Construct_opposite_direction_2;

  typedef Construct_direction_at_endpoint_2<Kernel, X_curve>
                                              Construct_direction_2;
  typedef Compare_y_at_x_for_segments_2<Kernel, X_curve>
                                              Compare_y_at_x_2;

  inline Construct_direction_2 construct_direction_2_object() const
  { return Construct_direction_2(); }

  inline Compare_y_at_x_2 compare_y_at_x_2_object() const
  { return Compare_y_at_x_2(); }

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
  
  // Compare the y value of two curves in an epsilon environment to
  // the left of the x value of the input point
  Comparison_result 
  curve_compare_at_x_left(const X_curve &cv1, const X_curve &cv2, 
                          const Point_2 &q) const 
  {
    // If one of the curves is vertical then return EQUAL.
    if ( is_vertical_2_object()(cv1) || (is_vertical_2_object()(cv2)) ) 
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
    if ( is_vertical_2_object()(cv1) || (is_vertical_2_object()(cv2))  ) 
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
  

protected:
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
  
};

CGAL_END_NAMESPACE

#endif // CGAL_PM_SEGMENT_EXACT_TRAITS_H
