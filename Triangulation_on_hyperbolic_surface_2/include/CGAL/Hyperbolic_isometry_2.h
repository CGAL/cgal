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
// Author(s)     : Vincent Despré, Loïc Dubois, Marc Pouget, Monique Teillaud

#ifndef CGAL_HYPERBOLIC_ISOMETRY_2_H
#define CGAL_HYPERBOLIC_ISOMETRY_2_H

#include <CGAL/license/Triangulation_on_hyperbolic_surface_2.h>

#include <CGAL/Complex_number.h>

namespace CGAL {

/*
Represents a hyperbolic isometry in the Poincare disk model the hyperbolic plane.
The isometry f is stored as list (c0, c1, c2, c3) of 4 complex numbers,
so that f(z) = (c0 z + c1) / (c2 z + c3) holds on every complex z in the open unit disk.
*/
template<class Traits>
class Hyperbolic_isometry_2
{
public:
  typedef Hyperbolic_isometry_2<Traits>                    Self;
  typedef typename Traits::FT                              FT;
  typedef typename Traits::Complex                         Complex_number;
  typedef typename Traits::Hyperbolic_point_2              Point;

  Hyperbolic_isometry_2();
  Hyperbolic_isometry_2(const Complex_number& c0, const Complex_number& c1, const Complex_number& c2, const Complex_number& c3);

  void set_to_identity();

  // Can be used to set the coefficients manually. Warning: the implementation does not check that the resulting moebius transform fixes the unit circle.
  void set_coefficients(const Complex_number& c0, const Complex_number& c1, const Complex_number& c2, const Complex_number& c3);
  void set_coefficient(int index, const Complex_number& coefficient);

  // returns the index-th coefficient
  const Complex_number& get_coefficient(int index) const;

  // evaluates the isometry at point
  Point evaluate(const Point& point) const;
  Point operator()(const Point& point) const;

private:
  Complex_number coefficients_[4];
};

// returns the composition of two isometries.
template<class Traits>
  Hyperbolic_isometry_2<Traits>  operator*(const Hyperbolic_isometry_2<Traits>& iso1, const Hyperbolic_isometry_2<Traits>& iso2);

// template<class Traits> std::ostream& operator<<(std::ostream& s, const Hyperbolic_isometry_2<Traits>& isometry);

// When inverse is 'false', returns the hyperbolic translation that maps -p to zero, and zero to p.
// Otherwise, returns the hyperbolic translation that maps p to zero, and zero to -p.
template<class Traits>
Hyperbolic_isometry_2<Traits> hyperbolic_translation(const typename Traits::Hyperbolic_point_2& p,
                                                     bool inverse = false);

// When inverse is 'false', returns the hyperbolic rotation around zero that maps p to q.
// Otherwise, returns the hyperbolic rotation around zero that maps q to p.
template<class Traits>
Hyperbolic_isometry_2<Traits> hyperbolic_rotation(const typename Traits::Hyperbolic_point_2& p,
                                                  const typename Traits::Hyperbolic_point_2& q,
                                                  bool inverse = false);

// returns the hyperbolic isometry that maps p1 to q1 and p2 to q2
template<class Traits>
Hyperbolic_isometry_2<Traits> isometry_pairing_the_sides(const typename Traits::Hyperbolic_point_2& p1,
                                                         const typename Traits::Hyperbolic_point_2& p2,
                                                         const typename Traits::Hyperbolic_point_2& q1,
                                                         const typename Traits::Hyperbolic_point_2& q2);

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template<class Traits>
Hyperbolic_isometry_2<Traits>::
Hyperbolic_isometry_2()
{
  set_to_identity();
}

template<class Traits>
Hyperbolic_isometry_2<Traits>::
Hyperbolic_isometry_2(const Complex_number& c0,
                      const Complex_number& c1,
                      const Complex_number& c2,
                      const Complex_number& c3)
{
  set_coefficients(c0,c1,c2,c3);
}
////////////////////////////////////////////////////////////////////////////////

template<class Traits>
void
Hyperbolic_isometry_2<Traits>::
set_to_identity()
{
  set_coefficients(Complex_number(FT(1)),
                   Complex_number(FT(0)),
                   Complex_number(FT(0)),
                   Complex_number(FT(1)));
}

template<class Traits>
void
Hyperbolic_isometry_2<Traits>::
set_coefficients(const Complex_number& c0,
                 const Complex_number& c1,
                 const Complex_number& c2,
                 const Complex_number& c3)
{
  set_coefficient(0, c0);
  set_coefficient(1, c1);
  set_coefficient(2, c2);
  set_coefficient(3, c3);
}

template<class Traits>
void
Hyperbolic_isometry_2<Traits>::
set_coefficient(int index, const Complex_number& coefficient)
{
  coefficients_[index] = coefficient;
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits>
const typename Traits::Complex&
Hyperbolic_isometry_2<Traits>::
get_coefficient(int index) const
{
  return coefficients_[index];
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits>
typename Traits::Hyperbolic_point_2
Hyperbolic_isometry_2<Traits>::
evaluate(const Point& point) const
{
  Complex_number z(point.x(), point.y());
  Complex_number numerator_of_the_result = coefficients_[0] * z + coefficients_[1];
  Complex_number denominator_of_the_result = coefficients_[2] * z + coefficients_[3];
  Complex_number result = numerator_of_the_result / denominator_of_the_result;
  return Point(result.real(), result.imag());
}

template<class Traits>
typename Traits::Hyperbolic_point_2
Hyperbolic_isometry_2<Traits>::
operator()(const Point& point) const
{
  return evaluate(point);
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits>
Hyperbolic_isometry_2<Traits> operator*(const Hyperbolic_isometry_2<Traits>& iso1,
                                        const Hyperbolic_isometry_2<Traits>& iso2)
{
  Hyperbolic_isometry_2<Traits> result;
  for (int i=0; i<2; i++) {
    for (int j=0; j<2; j++) {
      result.set_coefficient(2*i+j,
        iso1.get_coefficient(2*i) * iso2.get_coefficient(j) + iso1.get_coefficient(2*i+1) * iso2.get_coefficient(2+j));
    }
  }
  return result;
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits>
Hyperbolic_isometry_2<Traits> hyperbolic_translation(const typename Traits::Hyperbolic_point_2& p,
                                                     bool inverse)
{
  typename Traits::Complex one (typename Traits::FT(1));
  typename Traits::Complex z;
  if (inverse) {
    z = typename Traits::Complex(p.x(), p.y());
  } else {
    z = - typename Traits::Complex(p.x(), p.y());
  }
  Hyperbolic_isometry_2<Traits> result;
  result.set_coefficients(one, z, conj(z), one);
  return result;
}

template<class Traits>
Hyperbolic_isometry_2<Traits> hyperbolic_rotation(const typename Traits::Hyperbolic_point_2& p,
                                                  const typename Traits::Hyperbolic_point_2& q,
                                                  bool inverse)
{
  typename Traits::Complex zero (typename Traits::FT(0));

  Hyperbolic_isometry_2<Traits> result;
  if (inverse) {
    result.set_coefficients(typename Traits::Complex(p.x(), p.y()), zero, zero, typename Traits::Complex(q.x(), q.y()));
  } else {
    result.set_coefficients(typename Traits::Complex(q.x(), q.y()), zero, zero, typename Traits::Complex(p.x(), p.y()));
  }
  return result;
}

template<class Traits>
Hyperbolic_isometry_2<Traits> isometry_pairing_the_sides(const typename Traits::Hyperbolic_point_2& p1,
                                                         const typename Traits::Hyperbolic_point_2& p2,
                                                         const typename Traits::Hyperbolic_point_2& q1,
                                                         const typename Traits::Hyperbolic_point_2& q2)
{
  Hyperbolic_isometry_2<Traits> A,B,Binv,C;
  A = hyperbolic_translation<Traits>(p1);
  B = hyperbolic_translation<Traits>(q1);
  Binv = hyperbolic_translation<Traits>(q1, true);
  C = hyperbolic_rotation<Traits>(A.evaluate(p2), B.evaluate(q2));
  return (Binv*C)*A;
}

} // namespace CGAL

#endif // CGAL_HYPERBOLIC_ISOMETRY_2_H
