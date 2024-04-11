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

// This file contains the declaration and the implementation of the class Hyperbolic_surfaces_traits_2

#ifndef CGAL_HYPERBOLIC_SURFACES_TRAITS_2
#define CGAL_HYPERBOLIC_SURFACES_TRAITS_2

#include "Complex_without_sqrt.h"
#include <iostream>

namespace CGAL {

/*

template<class ComplexType>
class Hyperbolic_point_2 {
private:
  ComplexType _z;

public:
  typedef ComplexType                           Complex;

  Hyperbolic_point_2() { _z = Complex(); }
  Hyperbolic_point_2(const Complex& z) { _z = z; }

  Complex get_z() const { return _z; }
  void set_z(const Complex& z) { _z = z; }
};

template<class ComplexType>
std::ostream& operator<<(std::ostream& s, const Hyperbolic_point_2<ComplexType>& point) {
  s << point.get_z();
  return s;
}

template<class ComplexType>
void operator>>(std::istream& s, Hyperbolic_point_2<ComplexType>& point) {
  ComplexType z;
  s >> z;
  point.set_z(z);
}

*/


/*
Traits class offered by CGAL.
*/

/*
template<class FieldType>
class Hyperbolic_surfaces_traits_2 {
public:
  typedef FieldType                          FT;
  typedef Complex_without_sqrt<FieldType>    Complex;
  typedef Hyperbolic_point_2<Complex>        Point_2;
};
*/

template<class HyperbolicTraitsClass>
class Hyperbolic_surfaces_traits_2 : public HyperbolicTraitsClass {
public:
  typedef typename HyperbolicTraitsClass::FT                          FT;
  typedef typename HyperbolicTraitsClass::Hyperbolic_point_2          Hyperbolic_point_2;
  typedef Complex_without_sqrt<FT>                                    Complex;
};

} // namespace CGAL

#endif // CGAL_HYPERBOLIC_SURFACES_TRAITS_2
