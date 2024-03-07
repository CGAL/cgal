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

#include <CGAL/Complex_without_sqrt.h>
#include <CGAL/Hyperbolic_fundamental_domain_2.h>
#include <CGAL/Hyperbolic_surfaces_traits_2.h>
#include <iostream>

#include <CGAL/Gmpq.h>

using namespace CGAL;

typedef Gmpq FT;

typedef Hyperbolic_surfaces_traits_2<FT>                            Traits;
typedef Hyperbolic_fundamental_domain_2<Traits>                     Domain;

typedef typename Traits::Point_2                                    Point;
typedef Complex_without_sqrt<FT>                                    Complex;


int main() {
  std::vector<Point> vertices;
  Complex z0 = Complex (FT(809,10000),FT(0));
  Complex z1 = Complex (FT(7359,10000),FT(1877,10000));
  Complex z2 = Complex (FT(-999,2500),FT(881,1000));
  Complex z3 = Complex (FT("-22088524601252853411192791001942853611410938513/24711029456888649611435724068315791591836010000"),FT("9482675065452890527617859332378101016513362487/24711029456888649611435724068315791591836010000"));
  vertices.push_back( Point(z0) );
  vertices.push_back( Point(z1) );
  vertices.push_back( Point(z2) );
  vertices.push_back( Point(z3) );
  vertices.push_back( Point(-z0) );
  vertices.push_back( Point(-z1) );
  vertices.push_back( Point(-z2) );
  vertices.push_back( Point(-z3) );

  std::vector<int> pairings;
  for (int k=0; k<8; k++){
    pairings.push_back((k+4)%8);
  }

  Domain domain = Domain(vertices, pairings);
  assert( domain.size()==8 );
  for (int k=0; k<8; k++){
    assert( domain.vertex(k).get_z()==vertices[k].get_z());
    assert( domain.paired_side(k)==(k+4)%8 );

    assert( domain.side_pairing(k).evaluate(domain.vertex((k+4)%8)).get_z()==domain.vertex((k+1)%8).get_z() );
    assert( domain.side_pairing(k).evaluate(domain.vertex((k+5)%8)).get_z()==domain.vertex(k).get_z() );
  }

  assert( domain.is_valid() );

  Domain domain_prime;
  domain_prime.set(vertices, pairings);
  assert( domain_prime.size()==8 );
  for (int k=0; k<8; k++){
    assert( domain_prime.vertex(k).get_z()==vertices[k].get_z());
    assert( domain_prime.paired_side(k)==(k+4)%8 );

    assert( domain_prime.side_pairing(k).evaluate(domain_prime.vertex((k+4)%8)).get_z()==domain_prime.vertex((k+1)%8).get_z() );
    assert( domain_prime.side_pairing(k).evaluate(domain_prime.vertex((k+5)%8)).get_z()==domain_prime.vertex(k).get_z() );
  }

  Domain domain_ter = Domain();
  std::stringstream buffer;
  buffer << domain;
  buffer >> domain_ter;
  assert( domain_ter.size()==8 );
  for (int k=0; k<8; k++){
    assert( domain_ter.vertex(k).get_z()==vertices[k].get_z());
    assert( domain_ter.paired_side(k)==(k+4)%8 );
  }

  std::cout << "printing a domain for test purposes : " << std::endl << domain << std::endl;

  return 0;
}
