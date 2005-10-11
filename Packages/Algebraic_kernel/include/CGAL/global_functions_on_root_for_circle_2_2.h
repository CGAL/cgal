#ifndef CGAL_ALGEBRAIC_KERNEL_GLOBAL_FUNCTIONS_ON_ROOT_FOR_CIRCLE_2_2_H
#define CGAL_ALGEBRAIC_KERNEL_GLOBAL_FUNCTIONS_ON_ROOT_FOR_CIRCLE_2_2_H

CGAL_BEGIN_NAMESPACE

template < class AK >
inline 
Comparison_result 
compare_x(const typename AK::Root_for_circles_2_2& r1,
	   const typename AK::Root_for_circles_2_2& r2)
{ return AK().compare_x_object()(r1, r2); }

template < class AK >
inline 
Comparison_result 
compare_y(const typename AK::Root_for_circles_2_2& r1,
	   const typename AK::Root_for_circles_2_2& r2)
{ return AK().compare_y_object()(r1, r2); }

template < class AK >
inline 
Comparison_result 
compare_xy(const typename AK::Root_for_circles_2_2& r1,
	     const typename AK::Root_for_circles_2_2& r2)
{ return AK().compare_xy_object()(r1, r2); }

CGAL_END_NAMESPACE

#endif //CGAL_ALGEBRAIC_KERNEL_GLOBAL_FUNCTIONS_ON_ROOT_FOR_CIRCLE_2_2_H
