#ifndef CGAL_ALGEBRAIC_KERNEL_FUNCTION_OBJECTS_ON_ROOTS_AND_POLYNOMIALS_2_H
#define CGAL_ALGEBRAIC_KERNEL_FUNCTION_OBJECTS_ON_ROOTS_AND_POLYNOMIALS_2_H

#include <CGAL/Algebraic_kernel/internal_functions_on_roots_and_polynomials_2_2.h>

namespace CGAL {
  namespace AlgebraicFunctors {

  template < class AK >
  class Solve
  {
    typedef typename AK::Polynomial_for_circles_2_2 Equation;

    public:

    template < class OutputIterator >
    OutputIterator
    operator()(const Equation & e1, const Equation & e2,OutputIterator res)
      { return solve<AK> ( e1, e2, res); }

  };

  } // namespace AlgebraicFunctors
} // namespace CGAL

#endif // CGAL_ALGEBRAIC_KERNEL_FUNCTION_OBJECTS_ON_ROOTS_AND_POLYNOMIALS_2_H
