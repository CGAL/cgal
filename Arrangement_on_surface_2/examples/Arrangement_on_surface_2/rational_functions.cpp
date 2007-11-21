//! \file examples/Arrangement_2/ex_rational_functions.cpp
// Constructing an arrangement of arcs of rational functions.
#include <CGAL/basic.h>

#ifndef CGAL_USE_CORE
#include <iostream>
int main ()
{
  std::cout << "Sorry, this example needs CORE ..." << std::endl;
  return (0);
}
#else

#include <CGAL/Cartesian.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_rational_arc_traits_2.h>
#include <CGAL/Arrangement_2.h>

typedef CGAL::CORE_algebraic_number_traits            Nt_traits;
typedef Nt_traits::Rational                           Rational;
typedef Nt_traits::Algebraic                          Algebraic;
typedef CGAL::Cartesian<Algebraic>                    Alg_kernel;
typedef CGAL::Arr_rational_arc_traits_2<Alg_kernel,
                                        Nt_traits>    Traits_2;
typedef Traits_2::Point_2                             Point_2;
typedef Traits_2::Curve_2                             Rational_arc_2;
typedef Traits_2::Rat_vector                          Rat_vector;
typedef std::list<Rational_arc_2>                     Rat_arcs_list;
typedef CGAL::Arrangement_2<Traits_2>                 Arrangement_2;

int main ()
{
  // Create an arc supported by the polynomial y = x^4 - 6x^2 + 8,
  // defined over the interval [-2.1, 2.1]:
  Rat_vector        P1(5);
  P1[4] = 1; P1[3] = 0; P1[2] = -6; P1[1] = 0; P1[0] = 8;

  Rational_arc_2    a1 (P1, Algebraic(-2.1), Algebraic(2.1));

  // Create an arc supported by the function y = x / (1 + x^2),
  // defined over the interval [-3, 3]:
  Rat_vector        P2(2);
  P2[1] = 1; P2[0] = 0;

  Rat_vector        Q2(3);
  Q2[2] = 1; Q2[1] = 0; Q2[0] = 1;

  Rational_arc_2    a2 (P2, Q2, Algebraic(-3), Algebraic(3));

  // Create an arc supported by the parbola y = 8 - x^2,
  // defined over the interval [-2, 3]:
  Rat_vector        P3(5);
  P3[2] = -1; P3[1] = 0; P3[0] = 8;

  Rational_arc_2    a3 (P3, Algebraic(-2), Algebraic(3));

  // Create an arc supported by the line y = -2x,
  // defined over the interval [-3, 0]:
  Rat_vector        P4(2);
  P4[1] = -2; P4[0] = 0;

  Rational_arc_2    a4 (P4, Algebraic(-3), Algebraic(0));

  // Construct the arrangement of the four arcs.
  Arrangement_2              arr;
  std::list<Rational_arc_2>  arcs;

  arcs.push_back (a1);
  arcs.push_back (a2);
  arcs.push_back (a3);
  arcs.push_back (a4);
  insert (arr, arcs.begin(), arcs.end());

  // Print the arrangement size.
  std::cout << "The arrangement size:" << std::endl
            << "   V = " << arr.number_of_vertices()
            << ",  E = " << arr.number_of_edges()
            << ",  F = " << arr.number_of_faces() << std::endl;

  return 0;
}

#endif
