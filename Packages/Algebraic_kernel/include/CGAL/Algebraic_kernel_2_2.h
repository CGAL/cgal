#ifndef CGAL_ALGEBRAIC_KERNEL_2_2_H
#define CGAL_ALGEBRAIC_KERNEL_2_2_H

#include <CGAL/Root_of_2.h>
#include <CGAL/Polynomials_2_2.h>
#include <CGAL/Polynomials_1_2.h>
#include <CGAL/Root_for_circles_2_2.h>

#include <CGAL/Algebraic_kernel/function_objects_on_roots_and_polynomials_2_2.h> 
#include <CGAL/global_functions_on_roots_and_polynomials_2_2.h>
#include <CGAL/global_functions_on_roots_and_polynomial_1_2_and_2_2.h>

namespace CGAL {

  template< class RT_ >
  struct Algebraic_kernel_2_2
  {
    typedef Algebraic_kernel_2_2<RT_>        Self;

    typedef RT_  RT;

    typedef typename Root_of_traits< RT >::RootOf_2         Root_of_2;
    typedef typename Root_of_traits< RT >::RootOf_1         FT;
    typedef Root_for_circles_2_2< RT >             Root_for_circles_2_2;

    typedef Polynomial_for_circles_2_2<RT>
    Polynomial_for_circles_2_2; // probleme RT / FT

    typedef Polynomial_1_2<RT>
    Polynomial_1_2; // probleme RT / FT

    typedef AlgebraicFunctors::Solve<Self>  Solve;
    typedef AlgebraicFunctors::Sign_at<Self> Sign_at;

    typedef AlgebraicFunctors::Construct_polynomial_circle_2_2<Self>
    Construct_polynomial_circle_2_2;

    Solve 
    solve_object() const
      { return Solve(); }

    Construct_polynomial_circle_2_2
    construct_polynomial_circle_2_2_object() const
    {
      return Construct_polynomial_circle_2_2();
    }

    Sign_at
    sign_at_object() const
    { return Sign_at();}

  };

} // namespace CGAL

#endif // CGAL_ALGEBRAIC_KERNEL_2_2_H
