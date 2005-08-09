#ifndef CGAL_CURVED_KERNEL_FUNCTIONS_ON_CIRCLE_2_H
#define CGAL_CURVED_KERNEL_FUNCTIONS_ON_CIRCLE_2_H

namespace CGAL {

  // Should we have an iterator based interface, or both ?
  template <class CK>
  typename CK::Circular_arc_endpoint_2
  x_critical_points(const Circle_2<CK> & c, bool i)
  {
  	return CircularFunctors::x_critical_points<CK>(c,i);
  }

  template <class CK>
  typename CK::Circular_arc_endpoint_2
  y_critical_points(const Circle_2<CK> & c, bool i)
  {
  	return CircularFunctors::y_critical_points<CK>(c,i);
  }

} // namespace CGAL

#endif // CGAL_CURVED_KERNEL_FUNCTIONS_ON_CIRCLE_2_H
