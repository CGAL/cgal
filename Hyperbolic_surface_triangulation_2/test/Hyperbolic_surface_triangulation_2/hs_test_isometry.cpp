// Copyright (c) 2024
// INRIA Nancy (France), and Université Gustave Eiffel Marne-la-Vallee (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later EOR LicenseRef-Commercial
//
// Author(s)     : Vincent Despré, Loïc Dubois, Monique Teillaud

#include <iostream>

#include <CGAL/Gmpq.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_traits_2.h>
#include <CGAL/Hyperbolic_surfaces_traits_2.h>
#include <CGAL/Complex_without_sqrt.h>
#include <CGAL/Hyperbolic_isometry_2.h>

using namespace CGAL;

typedef Cartesian<Gmpq>                                                 Kernel;
typedef Hyperbolic_Delaunay_triangulation_traits_2<Kernel>              ParentTraits;
typedef Hyperbolic_surfaces_traits_2<ParentTraits>                      Traits;
typedef Hyperbolic_isometry_2<Traits>                                   Isometry;

typedef typename Traits::FT                                             FT;
typedef typename Traits::Hyperbolic_point_2                             Point;
typedef typename Traits::Complex                                        Complex;


int main() {
  Isometry identity_1 = Isometry ();
  assert( identity_1.get_coefficient(0)==Complex(FT(1)) );
  assert( identity_1.get_coefficient(1)==Complex(FT(0)) );
  assert( identity_1.get_coefficient(2)==Complex(FT(0)) );
  assert( identity_1.get_coefficient(3)==Complex(FT(1)) );

  Isometry identity_2;
  identity_2.set_to_identity();
  assert( identity_2.get_coefficient(0)==Complex(FT(1)) );
  assert( identity_2.get_coefficient(1)==Complex(FT(0)) );
  assert( identity_2.get_coefficient(2)==Complex(FT(0)) );
  assert( identity_2.get_coefficient(3)==Complex(FT(1)) );

  Isometry f;
  f.set_coefficient(0, Complex(FT(-12,17), FT(1,3)));
  f.set_coefficient(1, Complex(FT(56,7), FT(21,5)));
  f.set_coefficient(2, Complex(FT(56,7), FT(-21,5)));
  f.set_coefficient(3, Complex(FT(-12,17), FT(-1,3)));
  assert( f.get_coefficient(0)==Complex(FT(-12,17), FT(1,3)) );
  assert( f.get_coefficient(1)==Complex(FT(56,7), FT(21,5)) );
  assert( f.get_coefficient(2)==Complex(FT(56,7), FT(-21,5)) );
  assert( f.get_coefficient(3)==Complex(FT(-12,17), FT(-1,3)) );

  Isometry g;
  g.set_coefficients(Complex(FT(-12,17), FT(1,3)), Complex(FT(56,7), FT(21,5)), Complex(FT(56,7), FT(-21,5)), Complex(FT(-12,17), FT(-1,3)));
  assert( g.get_coefficient(0)==Complex(FT(-12,17), FT(1,3)) );
  assert( g.get_coefficient(1)==Complex(FT(56,7), FT(21,5)) );
  assert( g.get_coefficient(2)==Complex(FT(56,7), FT(-21,5)) );
  assert( g.get_coefficient(3)==Complex(FT(-12,17), FT(-1,3)) );

  Isometry h = f.compose(g);
  assert( h.get_coefficient(0) == Complex(FT(5333816,65025),FT(-8,17)) );
  assert( h.get_coefficient(1) == Complex(FT(-192,17),FT(-504,85)) );
  assert( h.get_coefficient(2) == Complex(FT(-192,17),FT(504,85)) );
  assert( h.get_coefficient(3) == Complex(FT(5333816,65025),FT(8,17)) );

  Point point (FT(3,11),FT(-1,73));
  Point image_point = h.evaluate(point);
  assert( image_point==Point(FT(9146011623056232,66567955527962869), FT(-12617302915955411,133135911055925738)) );

  std::cout << "printing an isometry for test purposes : " << std::endl << h;

  Isometry tau_1 = hyperbolic_translation<Traits>(point);
  Isometry tau_1_prime = hyperbolic_translation<Traits>(Point (FT(-3,11),FT(1,73)), true);
  Isometry tau_1_inv = hyperbolic_translation<Traits>(point, true);
  assert( tau_1.evaluate(image_point) == tau_1_prime.evaluate(image_point) );
  assert( tau_1.compose(tau_1_inv).evaluate(image_point) == image_point );

  Point p (FT(2,15),FT(0));
  Point q (FT(0),FT(17,93));
  Isometry rotation = hyperbolic_rotation<Traits>(p, q);
  Isometry rotation_prime = hyperbolic_rotation<Traits>(q, p, true);
  Isometry rotation_inv = hyperbolic_rotation<Traits>(p, q, true);
  assert( rotation.evaluate(image_point) == rotation_prime.evaluate(image_point) );
  assert( rotation.compose(rotation_inv).evaluate(image_point) == image_point );

  Point p_imag = rotation.evaluate(p);
  Point q_imag = rotation.evaluate(q);
  Isometry pairing = isometry_pairing_the_sides<Traits>(p, q, p_imag, q_imag);
  assert( pairing.evaluate(p) == p_imag );
  assert( pairing.evaluate(q) == q_imag );

  return 0;
}
