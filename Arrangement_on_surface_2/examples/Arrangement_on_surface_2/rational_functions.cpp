//! \file examples/Arrangement_2/rational_functions.cpp
// Constructing an arrangement of arcs of rational functions.

#include <CGAL/config.h>

#ifndef CGAL_USE_CORE
#include <iostream>
int main ()
{
  std::cout << "Sorry, this example needs CORE ..." << std::endl;
  return 0;
}

#else

#include <CGAL/CORE_BigInt.h>                      // NT
#include <CGAL/Algebraic_kernel_d_1.h>             // Algebraic Kernel
#include <CGAL/Arr_rational_function_traits_2.h>   // Traits
#include <CGAL/Arrangement_2.h>                    // Arrangement

typedef CORE::BigInt                               Number_type;
typedef CGAL::Algebraic_kernel_d_1<Number_type>           AK1;
typedef CGAL::Arr_rational_function_traits_2<AK1>  Traits_2;

typedef Traits_2::Polynomial_1                     Polynomial_1;
typedef Traits_2::Algebraic_real_1                 Alg_real_1;

typedef CGAL::Arrangement_2<Traits_2>              Arrangement_2;

int main ()
{
  CGAL::set_pretty_mode(std::cout);             // for nice printouts.

  // create a polynomial representing x .-)
  Polynomial_1 x = CGAL::shift(Polynomial_1(1),1);

  // Traits class object
  Traits_2 traits;
  Traits_2::Construct_x_monotone_curve_2 construct_arc
    = traits.construct_x_monotone_curve_2_object();

  // container storing all arcs
  std::vector<Traits_2::X_monotone_curve_2>  arcs;

  // Create an arc supported by the polynomial y = x^4 - 6x^2 + 8,
  // defined over the interval [-2.1, 2.1]:
  Polynomial_1 P1 = x*x*x*x - 6*x*x + 8;
  Alg_real_1 l(Traits_2::Algebraic_kernel_d_1::Bound(-2.1));
  Alg_real_1 r(Traits_2::Algebraic_kernel_d_1::Bound(2.1));
  arcs.push_back(construct_arc(P1, l, r));

  // Create an arc supported by the function y = x / (1 + x^2),
  // defined over the interval [-3, 3]:
  Polynomial_1 P2 = x;
  Polynomial_1 Q2 = 1+x*x;

  arcs.push_back(construct_arc(P2, Q2, Alg_real_1(-3), Alg_real_1(3)));

  // Create an arc supported by the parabola y = 8 - x^2,
  // defined over the interval [-2, 3]:
  Polynomial_1 P3 = 8 - x*x;
  arcs.push_back(construct_arc(P3, Alg_real_1(-2), Alg_real_1(3)));

  // Create an arc supported by the line y = -2x,
  // defined over the interval [-3, 0]:
  Polynomial_1 P4 = -2*x;
  arcs.push_back(construct_arc(P4, Alg_real_1(-3), Alg_real_1(0)));

  // Construct the arrangement of the four arcs.

  // Print the arcs.
  for (unsigned int i(0); i < arcs.size(); ++i)
    std::cout << arcs[i]<<std::endl;


  Arrangement_2 arr(&traits);
  insert(arr, arcs.begin(), arcs.end());

  // Print the arrangement size.
  std::cout << "The arrangement size:" << std::endl
            << "   V = " << arr.number_of_vertices()
            << ",  E = " << arr.number_of_edges()
            << ",  F = " << arr.number_of_faces() << std::endl;

  return 0;
}

#endif
