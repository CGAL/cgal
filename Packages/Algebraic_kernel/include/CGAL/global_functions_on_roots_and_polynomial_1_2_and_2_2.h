#ifndef CGAL_CURVED_KERNEL_GLOBAL_FUNCTIONS_ON_ROOTS_AND_POLYNOMIAL_1_2_AND_2_2_H
#define CGAL_CURVED_KERNEL_GLOBAL_FUNCTIONS_ON_ROOTS_AND_POLYNOMIAL_1_2_AND_2_2_H

namespace CGAL {

template< class AK, class OutputIterator>
inline
OutputIterator
solve( const typename AK::Polynomial_1_2 & e1,
       const typename AK::Polynomial_for_circles_2_2 & e2,
       OutputIterator res )
{
  return AK().solve_object()(e1,e2,res);
}

  template < class AK >
    inline 
    Sign sign_at( const typename AK::Polynomial_1_2 & equation,
		  const typename AK::Root_for_circles_2_2 r){
    return AK().sign_at_object()(equation, r);
 }





} // namespace CGAL

#endif // CGAL_CURVED_KERNEL_GLOBAL_FUNCTIONS_ON_ROOTS_AND_POLYNOMIAL_1_2_AND_2_2_H
