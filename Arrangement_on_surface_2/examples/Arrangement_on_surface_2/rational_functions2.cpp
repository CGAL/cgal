//! \file examples/Arrangement_2/ex_rational_functions.cpp
// Constructing an arrangement of arcs of rational functions.

#include <CGAL/basic.h>

#ifndef CGAL_USE_CORE
#include <iostream>
int main ()
{
  std::cout << "Sorry, this example needs CORE ..." << std::endl; 
  return 0;
}

#else

#include <CGAL/CORE_BigInt.h>                      // Integer
#include <CGAL/CORE_BigRat.h>                      // Rational
#include <CGAL/Algebraic_kernel_d_1.h>             // Algebraic Kernel
#include <CGAL/Arr_rational_function_traits_2.h>   // Traits
#include <CGAL/Arrangement_2.h>                    // Arrangement

typedef CORE::BigInt                               Integer;
typedef CORE::BigRat                               Rational;
typedef CGAL::Algebraic_kernel_d_1<Integer>	   AK1;
typedef CGAL::Arr_rational_function_traits_2<AK1>  Traits_2;

typedef std::vector<Rational>                      Rat_vec;
typedef Traits_2::Algebraic_real_1                 Alg_real_1;

typedef CGAL::Arrangement_2<Traits_2>              Arrangement_2;

int main ()
{
  CGAL::set_pretty_mode(std::cout);             // for nice printouts.

  // Traits class object 
  Traits_2 traits; 
  Traits_2::Construct_x_monotone_curve_2 construct_arc
    = traits.construct_x_monotone_curve_2_object(); 

  // container storing all arcs 
  std::vector<Traits_2::X_monotone_curve_2>  arcs;
  
  // Create an arc supported by the function y = 0.1x^4 - 0.6x^2 + 0.8 / 0.1,
  // defined over the interval [-2.1, 2.1]:
  Rat_vec P1,Q1;
  P1.push_back(Rational(8,10));
  P1.push_back(Rational(0));
  P1.push_back(Rational(-6,10));
  P1.push_back(Rational(0));
  P1.push_back(Rational(1,10));
  
  Q1.push_back(Rational(1,10));
  
  Alg_real_1 l(Traits_2::Algebraic_kernel_d_1::Bound(-2.1));
  Alg_real_1 r(Traits_2::Algebraic_kernel_d_1::Bound(2.1));
  arcs.push_back(construct_arc(P1.begin(), P1.end(), Q1.begin(), Q1.end(), l, r));

  // Create an arc supported by the function y = 0.1x / (0.1 + 0.1x^2),
  // defined over the interval [-3, 3]:
  Rat_vec P2,Q2;
  P2.push_back(Rational(0));
  P2.push_back(Rational(1,10));
  
  Q2.push_back(Rational(1,10));
  Q2.push_back(Rational(0));
  Q2.push_back(Rational(1,10));
  
  arcs.push_back(construct_arc(P2.begin(), P2.end(), Q2.begin(), Q2.end(),
                               Alg_real_1(-3), Alg_real_1(3)));
  
  // Create an arc supported by the parbola y = 0.8 - 0.1x^2 / 0.1,
  // defined over the interval [-2, 3]:
  Rat_vec P3,Q3;
  P3.push_back(Rational(4,5));
  P3.push_back(Rational(0));
  P3.push_back(Rational(-1,10));

  Q3.push_back(Rational(1,10));
  
  arcs.push_back(construct_arc(P3.begin(), P3.end(), Q3.begin(), Q3.end(),
                               Alg_real_1(-2), Alg_real_1(3)));
  
  // Create an arc supported by the line y = -0.2x / 0.1,
  // defined over the interval [-3, 0]:
  Rat_vec P4,Q4;
  P4.push_back(Rational(0));
  P4.push_back(Rational(-1,5));
  Q4.push_back(Rational(1,10));
  arcs.push_back(construct_arc(P4.begin(), P4.end(), Q4.begin(), Q4.end(),
                               Alg_real_1(-3), Alg_real_1(0)));

  // Print the arcs.
  for (unsigned int i(0); i < arcs.size(); ++i)
    std::cout << arcs[i]<<std::endl;

  // Construct the arrangement of the four arcs.
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
