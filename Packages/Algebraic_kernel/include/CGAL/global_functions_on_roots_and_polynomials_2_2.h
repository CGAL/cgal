#ifndef CGAL_CURVED_KERNEL_GLOBAL_FUNCTIONS_ON_ROOTS_AND_POLYNOMIALS_2_2_H
#define CGAL_CURVED_KERNEL_GLOBAL_FUNCTIONS_ON_ROOTS_AND_POLYNOMIALS_2_2_H

namespace CGAL {

template< class AK, class OutputIterator>
inline
OutputIterator
solve( const typename AK::Polynomial_for_circles_2_2 & e1,
       const typename AK::Polynomial_for_circles_2_2 & e2,
       OutputIterator res )
{
  return AK().solve_object()(e1,e2,res);
}

} // namespace CGAL

#endif // CGAL_CURVED_KERNEL_GLOBAL_FUNCTIONS_ON_ROOTS_AND_POLYNOMIALS_2_2_H
