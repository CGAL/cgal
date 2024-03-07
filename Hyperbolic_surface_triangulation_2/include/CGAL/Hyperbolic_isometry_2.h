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
template<class GeometricTraits_2>
class Hyperbolic_isometry_2{
private:
  typedef Hyperbolic_isometry_2<GeometricTraits_2>                    _Self;
  typedef typename GeometricTraits_2::FT                              _FT;
  typedef Complex_without_sqrt<_FT>                                   _Cmplx;
  typedef typename GeometricTraits_2::Point_2                         _Point;

  _Cmplx _coefficients[4];

public:
  typedef GeometricTraits_2 Geometric_traits_2;

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

template<class GeometricTraits_2> std::ostream& operator<<(std::ostream& s, const Hyperbolic_isometry_2<GeometricTraits_2>& isometry);

// When inverse=false, returns the hyperbolic translation that maps -p to zero, and zero to p. Otherwise, returns the hyperbolic translation that maps p to zero, and zero to -p.
template<class GeometricTraits_2>
Hyperbolic_isometry_2<GeometricTraits_2> hyperbolic_translation(const typename GeometricTraits_2::Point_2& p, bool inverse=false);

// When inverse=false, returns the hyperbolic rotation around zero that maps p to q. Otherwise, returns the hyperbolic rotation around zero that maps q to p.
template<class GeometricTraits_2>
Hyperbolic_isometry_2<GeometricTraits_2> hyperbolic_rotation(const typename GeometricTraits_2::Point_2& p, const typename GeometricTraits_2::Point_2& q, bool inverse=false);

// Returns the hyperbolic isometry that maps p1 to q1 and p2 to q2
template<class GeometricTraits_2>
Hyperbolic_isometry_2<GeometricTraits_2> segments_pairing(const typename GeometricTraits_2::Point_2& p1, const typename GeometricTraits_2::Point_2& p2, const typename GeometricTraits_2::Point_2& q1, const typename GeometricTraits_2::Point_2& q2);

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template<class GeometricTraits_2>
Hyperbolic_isometry_2<GeometricTraits_2>::Hyperbolic_isometry_2(){
  set_to_identity();
}

////////////////////////////////////////////////////////////////////////////////

template<class GeometricTraits_2>
void Hyperbolic_isometry_2<GeometricTraits_2>::set_to_identity(){
  set_coefficients(_Cmplx(_FT(1)),
                   _Cmplx(_FT(0)),
                   _Cmplx(_FT(0)),
                   _Cmplx(_FT(1)));
}

template<class GeometricTraits_2>
void Hyperbolic_isometry_2<GeometricTraits_2>::set_coefficients(const _Cmplx& c0, const _Cmplx& c1, const _Cmplx& c2, const _Cmplx& c3){
  set_coefficient(0, c0);
  set_coefficient(1, c1);
  set_coefficient(2, c2);
  set_coefficient(3, c3);
}

template<class GeometricTraits_2>
void Hyperbolic_isometry_2<GeometricTraits_2>::set_coefficient(int index, const _Cmplx& coefficient){
  _coefficients[index] = coefficient;
}

////////////////////////////////////////////////////////////////////////////////

template<class GeometricTraits_2>
Complex_without_sqrt<typename GeometricTraits_2::FT> Hyperbolic_isometry_2<GeometricTraits_2>::get_coefficient(int index) const{
  return _coefficients[index];
}

////////////////////////////////////////////////////////////////////////////////

template<class GeometricTraits_2>
typename GeometricTraits_2::Point_2 Hyperbolic_isometry_2<GeometricTraits_2>::evaluate(const _Point& point) const{
  _Cmplx z = point.get_z();
  _Cmplx numerator_of_the_result = _coefficients[0] * z + _coefficients[1];
  _Cmplx denominator_of_the_result = _coefficients[2] * z + _coefficients[3];

  return _Point(numerator_of_the_result / denominator_of_the_result);
}

////////////////////////////////////////////////////////////////////////////////

template<class GeometricTraits_2>
Hyperbolic_isometry_2<GeometricTraits_2> Hyperbolic_isometry_2<GeometricTraits_2>::compose(const Hyperbolic_isometry_2<GeometricTraits_2>& other) const{
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

template<class GeometricTraits_2>
std::ostream& operator<<(std::ostream& s, const Hyperbolic_isometry_2<GeometricTraits_2>& isometry){
  for (int k=0; k<4; k++){
    s << isometry.get_coefficient(k);
  }
  return s;
}

////////////////////////////////////////////////////////////////////////////////

template<class GeometricTraits_2>
Hyperbolic_isometry_2<GeometricTraits_2> hyperbolic_translation(const typename GeometricTraits_2::Point_2& p, bool inverse){
  Complex_without_sqrt<typename GeometricTraits_2::FT> one (typename GeometricTraits_2::FT(1));
  Complex_without_sqrt<typename GeometricTraits_2::FT> z;
  if (inverse){
    z = p.get_z();
  } else {
    z = - p.get_z();
  }
  Hyperbolic_isometry_2<GeometricTraits_2> result;
  result.set_coefficients(one, z, z.conjugate(), one);
  return result;
}

template<class GeometricTraits_2>
Hyperbolic_isometry_2<GeometricTraits_2> hyperbolic_rotation(const typename GeometricTraits_2::Point_2& p, const typename GeometricTraits_2::Point_2& q, bool inverse){
  Complex_without_sqrt<typename GeometricTraits_2::FT> zero (typename GeometricTraits_2::FT(0));
  Hyperbolic_isometry_2<GeometricTraits_2> result;
  if (inverse){
    result.set_coefficients(p.get_z(), zero, zero, q.get_z());
  } else {
    result.set_coefficients(q.get_z(), zero, zero, p.get_z());
  }
  return result;
}

template<class GeometricTraits_2>
Hyperbolic_isometry_2<GeometricTraits_2> segments_pairing(const typename GeometricTraits_2::Point_2& p1, const typename GeometricTraits_2::Point_2& p2, const typename GeometricTraits_2::Point_2& q1, const typename GeometricTraits_2::Point_2& q2){
    Hyperbolic_isometry_2<GeometricTraits_2> A,B,Binv,C;
    A = hyperbolic_translation<GeometricTraits_2>(p1);
    B = hyperbolic_translation<GeometricTraits_2>(q1);
    Binv = hyperbolic_translation<GeometricTraits_2>(q1,true);
    C = hyperbolic_rotation<GeometricTraits_2>(A.evaluate(p2), B.evaluate(q2));
    return (Binv.compose(C)).compose(A);
}

} // namespace CGAL

#endif // CGAL_HYPERBOLIC_ISOMETRY_2
