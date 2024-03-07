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
#include <CGAL/Hyperbolic_fundamental_domain_factory_2.h>
#include <CGAL/Hyperbolic_surfaces_traits_2.h>

#include <CGAL/Gmpq.h>

using namespace CGAL;

typedef Gmpq FT;

typedef Hyperbolic_surfaces_traits_2<FT>                            Traits;
typedef Hyperbolic_fundamental_domain_2<Traits>                     Domain;
typedef Hyperbolic_fundamental_domain_factory_2<Traits>             Factory;

typedef typename Traits::Point_2                                    Point;
typedef Complex_without_sqrt<FT>                                    Complex;


int main() {
  Factory factory (3459);
  Domain domain = factory.generate_domain_g2();

  std::vector<Point> vertices;
  Complex z0 = Complex (FT("4881/5000"),FT("0"));
  Complex z1 = Complex (FT("9211/10000"),FT("2733/10000"));
  Complex z2 = Complex (FT("1709/5000"),FT("7253/10000"));
  Complex z3 = Complex (FT("-427262704257582473474868322141310044732400799603/1267155016747148041260345910894159385550919570000"),FT("582571804584198065321856347012850217722442509611/1267155016747148041260345910894159385550919570000"));
  vertices.push_back( Point(z0) );
  vertices.push_back( Point(z1) );
  vertices.push_back( Point(z2) );
  vertices.push_back( Point(z3) );
  vertices.push_back( Point(-z0) );
  vertices.push_back( Point(-z1) );
  vertices.push_back( Point(-z2) );
  vertices.push_back( Point(-z3) );

  assert( domain.size()==8 );
  for (int k=0; k<8; k++){
    assert( domain.vertex(k).get_z()==vertices[k].get_z());
    assert( domain.paired_side(k)==(k+4)%8 );
  }

  return 0;
}
