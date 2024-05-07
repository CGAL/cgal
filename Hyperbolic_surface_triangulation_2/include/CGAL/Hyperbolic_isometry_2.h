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

// This file contains the declaration and the implementation of the class Hyperbolic_isometry_2

#ifndef CGAL_HYPERBOLIC_ISOMETRY_2
#define CGAL_HYPERBOLIC_ISOMETRY_2

#include "Complex_without_sqrt.h"

namespace CGAL {

/*
Represents a hyperbolic isometry in the Poincare disk model the hyperbolic plane. The isometry f is stored as list (c0, c1, c2, c3) of 4 complex numbers,
so that f(z) = (c0 z + c1) / (c2 z + c3) holds on every complex z in the open unit disk.
*/
template<class Traits>
class Hyperbolic_isometry_2{
public:
  typedef Hyperbolic_isometry_2<Traits>                    Self;
  typedef typename Traits::FT                              FT;
  typedef typename Traits::Complex                         ComplexNumber;
  typedef typename Traits::Hyperbolic_point_2              Point;

  Hyperbolic_isometry_2();

  void set_to_identity();

  // Can be used to set the coefficients manually. Be careful when doing so : the implementation does not check that the resulting moebius transform fixes the unit circle.
  void set_coefficients(const ComplexNumber& c0, const ComplexNumber& c1, const ComplexNumber& c2, const ComplexNumber& c3);
  void set_coefficient(int index, const ComplexNumber& coefficient);

  // Returns the index-th coefficient
  ComplexNumber get_coefficient(int index) const;

  // Evaluates the isometry at point
  Point evaluate(const Point& point) const;

  // Returns the composition of self and other
  Self compose(const Self& other) const;

private:
  ComplexNumber _coefficients[4];
};

template<class Traits> std::ostream& operator<<(std::ostream& s, const Hyperbolic_isometry_2<Traits>& isometry);

// When inverse=false, returns the hyperbolic translation that maps -p to zero, and zero to p. Otherwise, returns the hyperbolic translation that maps p to zero, and zero to -p.
template<class Traits>
Hyperbolic_isometry_2<Traits> hyperbolic_translation(const typename Traits::Hyperbolic_point_2& p, bool inverse=false);

// When inverse=false, returns the hyperbolic rotation around zero that maps p to q. Otherwise, returns the hyperbolic rotation around zero that maps q to p.
template<class Traits>
Hyperbolic_isometry_2<Traits> hyperbolic_rotation(const typename Traits::Hyperbolic_point_2& p, const typename Traits::Hyperbolic_point_2& q, bool inverse=false);

// Returns the hyperbolic isometry that maps p1 to q1 and p2 to q2
template<class Traits>
Hyperbolic_isometry_2<Traits> isometry_pairing_the_sides(const typename Traits::Hyperbolic_point_2& p1, const typename Traits::Hyperbolic_point_2& p2, const typename Traits::Hyperbolic_point_2& q1, const typename Traits::Hyperbolic_point_2& q2);

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template<class Traits>
Hyperbolic_isometry_2<Traits>::Hyperbolic_isometry_2(){
  set_to_identity();
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits>
void Hyperbolic_isometry_2<Traits>::set_to_identity(){
  set_coefficients(ComplexNumber(FT(1)),
                   ComplexNumber(FT(0)),
                   ComplexNumber(FT(0)),
                   ComplexNumber(FT(1)));
}

template<class Traits>
void Hyperbolic_isometry_2<Traits>::set_coefficients(const ComplexNumber& c0, const ComplexNumber& c1, const ComplexNumber& c2, const ComplexNumber& c3){
  set_coefficient(0, c0);
  set_coefficient(1, c1);
  set_coefficient(2, c2);
  set_coefficient(3, c3);
}

template<class Traits>
void Hyperbolic_isometry_2<Traits>::set_coefficient(int index, const ComplexNumber& coefficient){
  _coefficients[index] = coefficient;
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits>
typename Traits::Complex Hyperbolic_isometry_2<Traits>::get_coefficient(int index) const{
  return _coefficients[index];
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits>
typename Traits::Hyperbolic_point_2 Hyperbolic_isometry_2<Traits>::evaluate(const Point& point) const{
  ComplexNumber z (point.x(), point.y());
  ComplexNumber numerator_of_the_result = _coefficients[0] * z + _coefficients[1];
  ComplexNumber denominator_of_the_result = _coefficients[2] * z + _coefficients[3];
  ComplexNumber result = numerator_of_the_result / denominator_of_the_result;
  return Point(result.real_part(), result.imaginary_part());
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits>
Hyperbolic_isometry_2<Traits> Hyperbolic_isometry_2<Traits>::compose(const Hyperbolic_isometry_2<Traits>& other) const{
  Self result;
  for (int i=0; i<2; i++){
    for (int j=0; j<2; j++){
      result.set_coefficient(2*i+j,
            get_coefficient(2*i) * other.get_coefficient(j) + get_coefficient(2*i+1) * other.get_coefficient(2+j));
    }
  }
  return result;
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits>
std::ostream& operator<<(std::ostream& s, const Hyperbolic_isometry_2<Traits>& isometry){
  for (int k=0; k<4; k++){
    s << isometry.get_coefficient(k);
  }
  return s;
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits>
Hyperbolic_isometry_2<Traits> hyperbolic_translation(const typename Traits::Hyperbolic_point_2& p, bool inverse){
  typename Traits::Complex one (typename Traits::FT(1));
  typename Traits::Complex z;
  if (inverse){
    z = typename Traits::Complex(p.x(), p.y());
  } else {
    z = - typename Traits::Complex(p.x(), p.y());
  }
  Hyperbolic_isometry_2<Traits> result;
  result.set_coefficients(one, z, z.conjugate(), one);
  return result;
}

template<class Traits>
Hyperbolic_isometry_2<Traits> hyperbolic_rotation(const typename Traits::Hyperbolic_point_2& p, const typename Traits::Hyperbolic_point_2& q, bool inverse){
  typename Traits::Complex zero (typename Traits::FT(0));
  Hyperbolic_isometry_2<Traits> result;
  if (inverse){
    result.set_coefficients(typename Traits::Complex(p.x(), p.y()), zero, zero, typename Traits::Complex(q.x(), q.y()));
  } else {
    result.set_coefficients(typename Traits::Complex(q.x(), q.y()), zero, zero, typename Traits::Complex(p.x(), p.y()));
  }
  return result;
}

template<class Traits>
Hyperbolic_isometry_2<Traits> isometry_pairing_the_sides(const typename Traits::Hyperbolic_point_2& p1, const typename Traits::Hyperbolic_point_2& p2, const typename Traits::Hyperbolic_point_2& q1, const typename Traits::Hyperbolic_point_2& q2){
    Hyperbolic_isometry_2<Traits> A,B,Binv,C;
    A = hyperbolic_translation<Traits>(p1);
    B = hyperbolic_translation<Traits>(q1);
    Binv = hyperbolic_translation<Traits>(q1,true);
    C = hyperbolic_rotation<Traits>(A.evaluate(p2), B.evaluate(q2));
    return (Binv.compose(C)).compose(A);
}

} // namespace CGAL

#endif // CGAL_HYPERBOLIC_ISOMETRY_2
