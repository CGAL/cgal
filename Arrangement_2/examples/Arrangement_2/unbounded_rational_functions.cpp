//! \file examples/Arrangement_2/ex_unbounded_rational_functions.cpp
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
  std::list<Rational_arc_2>  arcs;

  // Create the rational functions (y = 1 / x), and (y = -1 / x).
  Rat_vector        P1(1);
  P1[0] = 1;

  Rat_vector        Q1(2);
  Q1[1] = 1; Q1[0] = 0;

  arcs.push_back (Rational_arc_2 (P1, Q1));

  P1[0] = -1;
  arcs.push_back (Rational_arc_2 (P1, Q1));

  // Create a bounded segments of the parabolas (y = -4*x^2 + 3) and
  // (y = 4*x^2 - 3), defined over [-sqrt(3)/2, sqrt(3)/2].
  const Algebraic   half_sqrt3 = CORE::sqrt(Algebraic(3)) / 2;
  Rat_vector        P2(3);
  P2[2] = -4; P2[1] = 0; P2[0] = 3;

  arcs.push_back (Rational_arc_2 (P2,
                                  -half_sqrt3, half_sqrt3));

  P2[2] = 4;  P2[0] = -3;
  arcs.push_back (Rational_arc_2 (P2,
                                  -half_sqrt3, half_sqrt3));

  // Create the rational function (y = 1 / 2*x) for x > 0, and the
  // rational function (y = -1 / 2*x) for x < 0.
  Rat_vector        P3(1);
  P3[0] = 1;

  Rat_vector        Q3(2);
  Q3[1] = 2; Q3[0] = 0;

  arcs.push_back (Rational_arc_2 (P3, Q3, Algebraic(0), true));

  P3[0] = -1;
  arcs.push_back (Rational_arc_2 (P3, Q3, Algebraic(0), false));

  // Construct the arrangement of the six arcs.
  Arrangement_2              arr;

  insert_curves (arr, arcs.begin(), arcs.end());

  // Print the arrangement size.
  std::cout << "The arrangement size:" << std::endl
            << "   V = " << arr.number_of_vertices()
            << " (plus " << arr.number_of_vertices_at_infinity()
            << " at infinity)"
            << ",  E = " << arr.number_of_edges()
            << ",  F = " << arr.number_of_faces()
            << " (" << arr.number_of_unbounded_faces() << " unbounded)"
            << std::endl << std::endl;

  return (0);
}

#endif
