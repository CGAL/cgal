// Copyright (c) 2024
// INRIA Nancy (France), and Université Gustave Eiffel Marne-la-Vallee (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Vincent Despré, Loïc Dubois, Monique Teillaud

#include <iostream>

#include <CGAL/Gmpq.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_traits_2.h>
#include <CGAL/Hyperbolic_surfaces_traits_2.h>
#include <CGAL/Complex_without_sqrt.h>
#include <CGAL/Hyperbolic_fundamental_domain_2.h>

using namespace CGAL;

typedef Cartesian<Gmpq>                                                 Kernel;
typedef Hyperbolic_Delaunay_triangulation_traits_2<Kernel>              ParentTraits;
typedef Hyperbolic_surfaces_traits_2<ParentTraits>                      Traits;
typedef Hyperbolic_fundamental_domain_2<Traits>                         Domain;

typedef typename Traits::FT                                             FT;
typedef typename Traits::Hyperbolic_point_2                             Point;
typedef typename Traits::Complex                                        Complex;


int main() {
  std::vector<Point> vertices;
  Point z0 = Point(FT("4881/5000"),FT("0"));
  Point z1 = Point(FT("9211/10000"),FT("2733/10000"));
  Point z2 = Point(FT("1709/5000"),FT("7253/10000"));
  Point z3 = Point(FT("-427262704257582473474868322141310044732400799603/1267155016747148041260345910894159385550919570000"),FT("582571804584198065321856347012850217722442509611/1267155016747148041260345910894159385550919570000"));
  Point z4 = Point(FT("-4881/5000"),FT("0"));
  Point z5 = Point(FT("-9211/10000"),FT("-2733/10000"));
  Point z6 = Point(FT("-1709/5000"),FT("-7253/10000"));
  Point z7 = Point(FT("427262704257582473474868322141310044732400799603/1267155016747148041260345910894159385550919570000"),FT("-582571804584198065321856347012850217722442509611/1267155016747148041260345910894159385550919570000"));
  vertices.push_back(z0);
  vertices.push_back(z1);
  vertices.push_back(z2);
  vertices.push_back(z3);
  vertices.push_back(z4);
  vertices.push_back(z5);
  vertices.push_back(z6);
  vertices.push_back(z7);

  std::vector<int> pairings;
  for (int k=0; k<8; k++){
    pairings.push_back((k+4)%8);
  }

  Domain domain = Domain(vertices, pairings);
  assert( domain.size()==8 );
  for (int k=0; k<8; k++){
    assert( domain.vertex(k)==vertices[k] );
    assert( domain.paired_side(k)==(k+4)%8 );

    assert( domain.side_pairing(k).evaluate(domain.vertex((k+4)%8))==domain.vertex((k+1)%8) );
    assert( domain.side_pairing(k).evaluate(domain.vertex((k+5)%8))==domain.vertex(k) );
  }

  assert( domain.is_valid() );

  Domain domain_prime;
  domain_prime.set(vertices, pairings);
  assert( domain_prime.size()==8 );
  for (int k=0; k<8; k++){
    assert( domain_prime.vertex(k)==vertices[k]);
    assert( domain_prime.paired_side(k)==(k+4)%8 );

    assert( domain_prime.side_pairing(k).evaluate(domain_prime.vertex((k+4)%8))==domain_prime.vertex((k+1)%8) );
    assert( domain_prime.side_pairing(k).evaluate(domain_prime.vertex((k+5)%8))==domain_prime.vertex(k) );
  }

  Domain domain_ter = Domain();
  std::stringstream buffer;
  buffer << domain;
  buffer >> domain_ter;
  assert( domain_ter.size()==8 );
  for (int k=0; k<8; k++){
    assert( domain_ter.vertex(k)==vertices[k]);
    assert( domain_ter.paired_side(k)==(k+4)%8 );
  }

  std::cout << "printing a domain for test purposes : " << std::endl << domain << std::endl;

  return 0;
}
