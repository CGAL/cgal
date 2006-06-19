//! \file examples/Arrangement_2/ex_unbounded_rational_functions.C
// Constructing an arrangement of unbounded portions of rational functions.
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
  // Create an arc that corresponds to the entire graph of y = 1/2*x^2.
  Rat_vector        P1(3);
  P1[2] = Rational (1, 2); P1[1] = 0; P1[0] = 0;

  Rational_arc_2    a1 (P1);

  // Create an arc that corresponds to graph of y = x + 1, for x = [0, +oo).
  Rat_vector        P2(2);
  P2[1] = 1; P2[0] = 1;

  Rational_arc_2    a2 (P2, Algebraic(0), true);  // Directed from 0 to right.

  // Create an arc that corresponds to graph of y = 1, for x = (-oo, 0].
  Rat_vector        P3(1);
  P3[0] = 1;

  Rational_arc_2    a3 (P3, Algebraic(0), false); // Directed from 0 to left.

  // Create the entire rational function y = 1 / (x^2 - 8x + 15).
  /*
  Rat_vector        P2(1);
  P2[0] = 1;

  Rat_vector        Q2(3);
  Q2[2] = 1; Q2[1] = -8; Q2[0] = 15;

  Rational_arc_2    a2 (P2, Q2);
  */

  // Construct the arrangement of the four arcs.
  Arrangement_2              arr;

  insert_curve (arr, a1);
  insert_curve (arr, a2);
  insert_curve (arr, a3);

  // Print the arrangement size.
  std::cout << "The arrangement size:" << std::endl
            << "   V = " << arr.number_of_vertices()
            << " (plus " << arr.number_of_vertices_at_infinity()
            << " at infinity)"
            << ",  E = " << arr.number_of_edges() 
            << ",  F = " << arr.number_of_faces() 
            << " (" << arr.number_of_unbounded_faces() << " unbounded)"
            << std::endl << std::endl;

  return 0;
}

#endif
