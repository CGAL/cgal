#ifndef CGAL_ALGEBRAIC_KERNEL_2_2_H
#define CGAL_ALGEBRAIC_KERNEL_2_2_H

#include <CGAL/Root_of_2.h>
#include <CGAL/Polynomials_2_2.h>
#include <CGAL/Polynomials_1_2.h>

#include <CGAL/Algebraic_kernel/function_objects_on_roots_and_polynomials_2_2.h> 
#include <CGAL/global_functions_on_roots_and_polynomials_2_2.h>

namespace CGAL {

  template< class RT >
  struct Algebraic_kernel_2_2
  {
    typedef Algebraic_kernel_2_2<RT>        Self;

    typedef typename Root_of_traits< RT >::RootOf_2         Root_of_2;
    typedef typename Root_of_traits< RT >::RootOf_1         FT;

    typedef Polynomial_for_circles_2_2<RT>
    Polynomial_for_circles_2_2; // probleme RT / FT
    typedef Polynomial_1_2<RT>
    Polynomial_1_2; // probleme RT / FT

    typedef AlgebraicFunctors::Solve<Self>  Solve;

    Solve 
    solve_object() const
      { return Solve(); }
  };

} // namespace CGAL

#endif // CGAL_ALGEBRAIC_KERNEL_2_2_H
