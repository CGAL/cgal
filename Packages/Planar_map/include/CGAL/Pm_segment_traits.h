#ifndef CGAL_PM_SEGMENT_TRAITS_H
#define CGAL_PM_SEGMENT_TRAITS_H

CGAL_BEGIN_NAMESPACE

template <class Kernel_>
class Pm_segment_traits
{
public:

  typedef Kernel_                    Kernel;

  // traits objects
  typedef typename Kernel::Point_2   Point_2;
  typedef Point_2                    Point;
  typedef typename Kernel::Vector_2  Vector_2;
  typedef Vector_2                   Vector;
  typedef typename Kernel::Segment_2 X_curve;

  // Things I get from the kernel
  typedef typename Kernel::Is_vertical_2      Is_vertical_2;
  typedef typename Kernel::Compare_y_at_x_2   Compare_y_at_x_2;
  typedef typename Kernel::Construct_opposite_segment_2 
                                              Construct_opposite_segment_2;
  typedef typename Kernel::Counterclockwise_in_between_2
                                              Counterclockwise_in_between_2;
  typedef typename Kernel::Equal_2            Equal_2;
  typedef typename Kernel::Has_on_2           Has_on_2;
  typedef typename Kernel::Compare_x_2        Compare_x_2;
  typedef typename Kernel::Compare_y_2        Compare_y_2;

  // Add to kernel?
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
  
  // done
  Point curve_source(const X_curve & cv) const 
  { 
    return m_kernel.construct_vertex_2_object()(cv, 0);
  }
  
  // done
  Point curve_target(const X_curve & cv) const 
  {
    return m_kernel.construct_vertex_2_object()(cv, 1);
  }
  
  // done
  bool curve_is_vertical(const X_curve & cv) const 
  {
    return m_kernel.is_vertical_2_object()(cv);
  }

  // Intorduce Is_in_x_range_2 / Is_in_x_closed_range_2 ?
  bool curve_is_in_x_range(const X_curve & cv, const Point & q) const
  { 
    return !( is_right(q, rightmost(cv.source(), cv.target())) ||
              is_left(q, leftmost(cv.source(), cv.target()))	 );
  }

  // Intorduce Is_in_y_range_2 / Intorduce Is_in_y_closed_range_2 ?
  bool curve_is_in_y_range(const X_curve &cv, const Point & q) const
  { 
    bool r = !( is_lower(q, lowest(cv.source(), cv.target())) ||
                is_higher(q, highest(cv.source(), cv.target())) );
    return r;
  }
  
  // Introduce Point_status_2 / Curve_point_status_2 ?
  Curve_point_status 
  curve_get_point_status(const X_curve &cv, const Point & p) const
  {
    if (!curve_is_in_x_range(cv, p))
      return CURVE_NOT_IN_RANGE;
	if (!curve_is_vertical(cv))
	{
		int res = compare_y(p, curve_calc_point(cv, p));
		if (res == SMALLER) return UNDER_CURVE;
		if (res == LARGER)	return ABOVE_CURVE;
		//if (res == EQUAL) 
		return ON_CURVE;
	}
	else
	{
		if (is_lower(p,lowest(curve_source(cv),curve_target(cv))))
			return UNDER_CURVE;
		if (is_higher(p,highest(curve_source(cv),curve_target(cv))))
			return ABOVE_CURVE;
		// if (curve_is_in_y_range(cv,p))
		return ON_CURVE;
	}
  }

  // done  
  Comparison_result 
  curve_compare_at_x(const X_curve &cv1, const X_curve &cv2, const Point &q) 
    const 
  {
    typename Kernel::Line_2 l1 = m_kernel.construct_line_2_object()(cv1);
    typename Kernel::Line_2 l2 = m_kernel.construct_line_2_object()(cv2);
    return m_kernel.compare_y_at_x_2_object()(q, l1, l2);
  }
  
  // Introduce Compare_x_left_of_2 ?
  Comparison_result 
  curve_compare_at_x_left(const X_curve &cv1, const X_curve &cv2, 
                          const Point &q) const 
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
    return compare_value(curve_derivative(cv2), curve_derivative(cv1));
  }
  
  // Introduce Compare_x_right_of_2 ?
  // (Constructing the opposite sounds too costly)
  Comparison_result 
  curve_compare_at_x_right(const X_curve &cv1, const X_curve &cv2, 
			   const Point & q) const 
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
    return compare_value(curve_derivative(cv1), curve_derivative(cv2));
  }
  
  // done
  X_curve curve_flip(const X_curve &cv) const
  {
    return m_kernel.construct_opposite_segment_2_object()(cv);
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

  // done
  Comparison_result compare_x(const Point &p1, const Point &p2) const
  { return m_kernel.compare_x_2_object()(p1, p2); }

  // done
  Comparison_result compare_y(const Point &p1, const Point &p2) const
  { return m_kernel.compare_y_2_object()(p1, p2); }

  // done
  bool curve_is_same(const X_curve &cv1, const X_curve &cv2) const
  {
    return m_kernel.equal_2_object()(cv1, cv2);
  }

  // done
  bool is_point_on_curve(const X_curve &cv, const Point& p) const //check
  {
    return m_kernel.has_on_2_object()(cv, p);
  }

private:
  bool is_left(const Point &p1, const Point &p2) const 
  { return m_kernel.less_x_2_object()(p1, p2); }
  bool is_right(const Point &p1, const Point &p2) const 
  { return m_kernel.less_x_2_object()(p2, p1); }
  bool is_same_x(const Point &p1, const Point &p2) const 
  { return m_kernel.equal_x_object()(p1, p2); }
  bool is_lower(const Point &p1, const Point &p2) const 
  { return m_kernel.less_y_2_object()(p1, p2); }
  bool is_higher(const Point &p1, const Point &p2) const 
  { return m_kernel.less_y_2_object()(p2, p1); }
  bool is_same_y(const Point &p1, const Point &p2) const 
  { return m_kernel.equal_y_object()(p1, p2); }
  bool is_same(const Point &p1, const Point &p2) const
  {
    return (compare_x(p1, p2) == EQUAL) &&
      (compare_y(p1, p2) == EQUAL);
  }
  const Point& leftmost(const Point &p1, const Point &p2) const
  { return (is_left(p1, p2) ? p1 : p2); }

  const Point& rightmost(const Point &p1, const Point &p2) const
  { return (is_right(p1, p2) ? p1 : p2); }
  
  const Point& lowest(const Point &p1, const Point &p2) const
  { return (is_lower(p1, p2) ? p1 : p2); }
  
  const Point& highest(const Point &p1, const Point &p2) const
  { return (is_higher(p1, p2) ? p1 : p2); }
  
  // Comment this one ! ##############
  Point curve_calc_point(const X_curve &cv, const Point & q) const
  {
    if ( ! curve_is_in_x_range(cv, q) )
      return cv.source();
    
    if (curve_is_vertical(cv))
      return cv.source();
    
    const Point & a = cv.source();
    const Point & b = cv.target();
    return Point ((b.hx() * a.hw() - a.hx() * b.hw()) * q.hx() * a.hw(),
                  (b.hx() * a.hw() - a.hx() * b.hw()) * q.hw() * a.hy() + 
                  (b.hy() * a.hw() - a.hy() * b.hw()) * 
                  (q.hx() * a.hw() - a.hx() * q.hw()),  
                  (b.hx() * a.hw() - a.hx() * b.hw()) * q.hw() * a.hw());
  }
  
  typename Kernel::FT curve_derivative(const X_curve &cv) const
  {
    CGAL_assertion(!curve_is_vertical(cv));
    
    return ( (cv.target()).y() - cv.source().y()) / 
      (cv.target().x() - cv.source().x());
  }
  
  Comparison_result compare_value(const typename Kernel::FT &v1, 
				  const typename Kernel::FT &v2) const
  {
    typename Kernel::FT       delta = v1 - v2;
    const typename Kernel::FT zero(0);
    if (delta == zero)
      return EQUAL;
    if (zero < delta)
      return LARGER;
    else
      return SMALLER;
  }
};

CGAL_END_NAMESPACE

#endif // CGAL_PM_SEGMENT_EXACT_TRAITS_H
