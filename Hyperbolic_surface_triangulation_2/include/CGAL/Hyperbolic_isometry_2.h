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
template<class Traits_2>
class Hyperbolic_isometry_2{
private:
  typedef Hyperbolic_isometry_2<Traits_2>                    _Self;
  typedef typename Traits_2::FT                              _FT;
  typedef typename Traits_2::Complex                         _Cmplx;
  typedef typename Traits_2::Hyperbolic_point_2              _Point;

  _Cmplx _coefficients[4];

public:
  Hyperbolic_isometry_2();

  void set_to_identity();

  // Can be used to set the coefficients manually. Be careful when doing so : the implementation does not check that the resulting moebius transform fixes the unit circle.
  void set_coefficients(const _Cmplx& c0, const _Cmplx& c1, const _Cmplx& c2, const _Cmplx& c3);
  void set_coefficient(int index, const _Cmplx& coefficient);

  // Returns the index-th coefficient
  _Cmplx get_coefficient(int index) const;

  // Evaluates the isometry at point
  _Point evaluate(const _Point& point) const;

  // Returns the composition of self and other
  _Self compose(const _Self& other) const;
};

template<class Traits_2> std::ostream& operator<<(std::ostream& s, const Hyperbolic_isometry_2<Traits_2>& isometry);

// When inverse=false, returns the hyperbolic translation that maps -p to zero, and zero to p. Otherwise, returns the hyperbolic translation that maps p to zero, and zero to -p.
template<class Traits_2>
Hyperbolic_isometry_2<Traits_2> hyperbolic_translation(const typename Traits_2::Hyperbolic_point_2& p, bool inverse=false);

// When inverse=false, returns the hyperbolic rotation around zero that maps p to q. Otherwise, returns the hyperbolic rotation around zero that maps q to p.
template<class Traits_2>
Hyperbolic_isometry_2<Traits_2> hyperbolic_rotation(const typename Traits_2::Hyperbolic_point_2& p, const typename Traits_2::Hyperbolic_point_2& q, bool inverse=false);

// Returns the hyperbolic isometry that maps p1 to q1 and p2 to q2
template<class Traits_2>
Hyperbolic_isometry_2<Traits_2> isometry_pairing_the_sides(const typename Traits_2::Hyperbolic_point_2& p1, const typename Traits_2::Hyperbolic_point_2& p2, const typename Traits_2::Hyperbolic_point_2& q1, const typename Traits_2::Hyperbolic_point_2& q2);

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template<class Traits_2>
Hyperbolic_isometry_2<Traits_2>::Hyperbolic_isometry_2(){
  set_to_identity();
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits_2>
void Hyperbolic_isometry_2<Traits_2>::set_to_identity(){
  set_coefficients(_Cmplx(_FT(1)),
                   _Cmplx(_FT(0)),
                   _Cmplx(_FT(0)),
                   _Cmplx(_FT(1)));
}

template<class Traits_2>
void Hyperbolic_isometry_2<Traits_2>::set_coefficients(const _Cmplx& c0, const _Cmplx& c1, const _Cmplx& c2, const _Cmplx& c3){
  set_coefficient(0, c0);
  set_coefficient(1, c1);
  set_coefficient(2, c2);
  set_coefficient(3, c3);
}

template<class Traits_2>
void Hyperbolic_isometry_2<Traits_2>::set_coefficient(int index, const _Cmplx& coefficient){
  _coefficients[index] = coefficient;
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits_2>
typename Traits_2::Complex Hyperbolic_isometry_2<Traits_2>::get_coefficient(int index) const{
  return _coefficients[index];
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits_2>
typename Traits_2::Hyperbolic_point_2 Hyperbolic_isometry_2<Traits_2>::evaluate(const _Point& point) const{
  _Cmplx z (point.x(), point.y());
  _Cmplx numerator_of_the_result = _coefficients[0] * z + _coefficients[1];
  _Cmplx denominator_of_the_result = _coefficients[2] * z + _coefficients[3];
  _Cmplx result = numerator_of_the_result / denominator_of_the_result;
  return _Point(result.real(), result.imag());
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits_2>
Hyperbolic_isometry_2<Traits_2> Hyperbolic_isometry_2<Traits_2>::compose(const Hyperbolic_isometry_2<Traits_2>& other) const{
  _Self result;
  for (int i=0; i<2; i++){
    for (int j=0; j<2; j++){
      result.set_coefficient(2*i+j,
            get_coefficient(2*i) * other.get_coefficient(j) + get_coefficient(2*i+1) * other.get_coefficient(2+j));
    }
  }
  return result;
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits_2>
std::ostream& operator<<(std::ostream& s, const Hyperbolic_isometry_2<Traits_2>& isometry){
  for (int k=0; k<4; k++){
    s << isometry.get_coefficient(k);
  }
  return s;
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits_2>
Hyperbolic_isometry_2<Traits_2> hyperbolic_translation(const typename Traits_2::Hyperbolic_point_2& p, bool inverse){
  typename Traits_2::Complex one (typename Traits_2::FT(1));
  typename Traits_2::Complex z;
  if (inverse){
    z = typename Traits_2::Complex(p.x(), p.y());
  } else {
    z = - typename Traits_2::Complex(p.x(), p.y());
  }
  Hyperbolic_isometry_2<Traits_2> result;
  result.set_coefficients(one, z, z.conjugate(), one);
  return result;
}

template<class Traits_2>
Hyperbolic_isometry_2<Traits_2> hyperbolic_rotation(const typename Traits_2::Hyperbolic_point_2& p, const typename Traits_2::Hyperbolic_point_2& q, bool inverse){
  typename Traits_2::Complex zero (typename Traits_2::FT(0));
  Hyperbolic_isometry_2<Traits_2> result;
  if (inverse){
    result.set_coefficients(typename Traits_2::Complex(p.x(), p.y()), zero, zero, typename Traits_2::Complex(q.x(), q.y()));
  } else {
    result.set_coefficients(typename Traits_2::Complex(q.x(), q.y()), zero, zero, typename Traits_2::Complex(p.x(), p.y()));
  }
  return result;
}

template<class Traits_2>
Hyperbolic_isometry_2<Traits_2> isometry_pairing_the_sides(const typename Traits_2::Hyperbolic_point_2& p1, const typename Traits_2::Hyperbolic_point_2& p2, const typename Traits_2::Hyperbolic_point_2& q1, const typename Traits_2::Hyperbolic_point_2& q2){
    Hyperbolic_isometry_2<Traits_2> A,B,Binv,C;
    A = hyperbolic_translation<Traits_2>(p1);
    B = hyperbolic_translation<Traits_2>(q1);
    Binv = hyperbolic_translation<Traits_2>(q1,true);
    C = hyperbolic_rotation<Traits_2>(A.evaluate(p2), B.evaluate(q2));
    return (Binv.compose(C)).compose(A);
}

} // namespace CGAL

#endif // CGAL_HYPERBOLIC_ISOMETRY_2
