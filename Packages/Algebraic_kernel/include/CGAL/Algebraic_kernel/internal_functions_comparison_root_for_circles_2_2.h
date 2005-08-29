// Author : Julien HAZEBROUCK

#ifndef CGAL_ALGEBRAIC_KERNEL_COMPARISON_ROOT_FOR_CIRCLES_2_2_H
#define CGAL_ALGEBRAIC_KERNEL_COMPARISON_ROOT_FOR_CIRCLES_2_2_H


namespace CGAL {
  namespace AlgebraicFunctors{
  
    template <typename RT>
      Comparison_result 
      compare_x(const CGAL::Root_for_circles_2_2<RT>& r1, const CGAL::Root_for_circles_2_2<RT>& r2){
      return compare(r1.x(), r2.x());
    }
  
    template <typename RT>
      Comparison_result 
      compare_y(const CGAL::Root_for_circles_2_2<RT>& r1, const CGAL::Root_for_circles_2_2<RT>& r2){
      return compare(r1.y(), r2.y());
    }
  
    template <typename RT>
      Comparison_result
      compare_xy(const CGAL::Root_for_circles_2_2<RT>& r1, const CGAL::Root_for_circles_2_2<RT>& r2){
      Comparison_result compx = compare_x(r1, r2);
      if(compx != 0)
	return compx;
      return compare_y(r1, r2);
    }

  } // namespace AlgebraicFunctors
} // namespace CGAL

#endif // CGAL_ALGEBRAIC_KERNEL_COMPARISON_ROOT_FOR_CIRCLES_2_2_H
