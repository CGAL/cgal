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

#ifndef CGAL_HYPERBOLIC_FUNDAMENTAL_DOMAIN_FACTORY_2_H
#define CGAL_HYPERBOLIC_FUNDAMENTAL_DOMAIN_FACTORY_2_H

#include <CGAL/license/Triangulation_on_hyperbolic_surface_2.h>

#include <CGAL/Hyperbolic_fundamental_domain_2.h>
#include <CGAL/Random.h>

#include <cmath>
#include <vector>

namespace CGAL {

/*
Factory class, whose only purpose is to construct random fundamental domains of
closed orientable hyperbolic surfaces. The function
`make_hyperbolic_fundamental_domain_g2()` constructs such a domain for a surface of
genus 2.
*/
template<class Traits>
class Hyperbolic_fundamental_domain_factory_2
{
private:
  typedef typename Traits::FT                    _FT;
  typedef typename Traits::Complex               _Cmplx;
  typedef typename Traits::Hyperbolic_point_2    _Point;

  Random _random;

public:
  Hyperbolic_fundamental_domain_factory_2();
  Hyperbolic_fundamental_domain_2<Traits> make_hyperbolic_fundamental_domain_g2(unsigned int seed);

private:
  float random_positive_float(); // returns number in [0,1]
  float random_float(); // returns number in [-1,1]
  Complex_number<float> random_complex_float(); // returns complex z such that modulus(z) < 1 and imag(z) > 0

  _FT exact_number_from_float(float x);
  _Cmplx exact_complex_from_float_complex(const Complex_number<float>& z);

  bool try_to_compute_inexact_z0_from_z1_z2_z3(Complex_number<float>& z0, Complex_number<float>& z1, Complex_number<float>& z2, Complex_number<float>& z3);
  bool try_to_compute_exact_z3_from_z0_z1_z2(_Cmplx& z0, _Cmplx& z1, _Cmplx& z2, _Cmplx& z3);

  bool sanity_check(_Cmplx& z0, _Cmplx& z1, _Cmplx& z2, _Cmplx& z3);

  const int CGAL_DENOMINATOR_FOR_GENERATION = 10000;
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template<class Traits>
  Hyperbolic_fundamental_domain_factory_2<Traits>::Hyperbolic_fundamental_domain_factory_2() {}

////////////////////////////////////////////////////////////////////////////////

template<class Traits>
Hyperbolic_fundamental_domain_2<Traits>
Hyperbolic_fundamental_domain_factory_2<Traits>::
make_hyperbolic_fundamental_domain_g2(unsigned int seed)
{
  _random = Random(seed);
  bool is_domain_generated = false;
  _Cmplx exact_z0, exact_z1, exact_z2, exact_z3;

  while (!is_domain_generated) {
    // 1. Generate inexact z0,z1,z2,z3
    Complex_number<float> z0, z1, z2, z3;
    z1 = random_complex_float();
    z2 = random_complex_float();
    z3 = random_complex_float();
    while (! try_to_compute_inexact_z0_from_z1_z2_z3(z0,z1,z2,z3)) {
      z1 = random_complex_float();
      z2 = random_complex_float();
      z3 = random_complex_float();
    }

    // 2. Compute exact z0,z1,z2,z3 nearby
    exact_z0 = exact_complex_from_float_complex(z0);
    exact_z1 = exact_complex_from_float_complex(z1);
    exact_z2 = exact_complex_from_float_complex(z2);
    exact_z3 = exact_complex_from_float_complex(z3);

    // 3. Modify z3 to fix the area...
    is_domain_generated = try_to_compute_exact_z3_from_z0_z1_z2(exact_z0, exact_z1, exact_z2, exact_z3);
    if (is_domain_generated) {
      // ... and perform a sanity check
      is_domain_generated = sanity_check(exact_z0, exact_z1, exact_z2, exact_z3);
    }
  }

  _Cmplx exact_zero(_FT(0), _FT(0));
  std::vector<_Point> vertices;
  vertices.push_back(_Point(exact_z0.real(), exact_z0.imag()));
  vertices.push_back(_Point(exact_z1.real(), exact_z1.imag()));
  vertices.push_back(_Point(exact_z2.real(), exact_z2.imag()));
  vertices.push_back(_Point(exact_z3.real(), exact_z3.imag()));
  vertices.push_back(_Point(-exact_z0.real(), -exact_z0.imag()));
  vertices.push_back(_Point(-exact_z1.real(), -exact_z1.imag()));
  vertices.push_back(_Point(-exact_z2.real(), -exact_z2.imag()));
  vertices.push_back(_Point(-exact_z3.real(), -exact_z3.imag()));

  std::vector<int> pairings;
  for (int k=0; k<8; ++k) {
    pairings.push_back((k+4)%8);
  }

  Hyperbolic_fundamental_domain_2<Traits> domain(vertices, pairings);
  return domain;
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits>
float Hyperbolic_fundamental_domain_factory_2<Traits>::random_positive_float() {
  return _random.uniform_01<float>();
}

template<class Traits>
float Hyperbolic_fundamental_domain_factory_2<Traits>::random_float() {
  return _random.uniform_01<float>() * 2  - 1;
}

template<class Traits>
Complex_number<float>
Hyperbolic_fundamental_domain_factory_2<Traits>::
random_complex_float()
{
  Complex_number<float> result (random_float(), random_positive_float());
  while (norm(result) >= 1) {
    result.real(random_float());
    result.imag(random_positive_float());
  }

  return result;
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits>
typename Traits::FT
Hyperbolic_fundamental_domain_factory_2<Traits>::
exact_number_from_float(float x)
{
  if (x < 0) {
    return _FT(0)-exact_number_from_float(-x);
  }
  return _FT(int(x*CGAL_DENOMINATOR_FOR_GENERATION)%CGAL_DENOMINATOR_FOR_GENERATION)/_FT(CGAL_DENOMINATOR_FOR_GENERATION);
}

template<class Traits>
typename Traits::Complex
Hyperbolic_fundamental_domain_factory_2<Traits>::
exact_complex_from_float_complex(const Complex_number<float>& z)
{
  return _Cmplx(exact_number_from_float(z.real()), exact_number_from_float(z.imag()));
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits>
bool Hyperbolic_fundamental_domain_factory_2<Traits>::
try_to_compute_inexact_z0_from_z1_z2_z3(Complex_number<float>& z0,
                                        Complex_number<float>& z1,
                                        Complex_number<float>& z2,
                                        Complex_number<float>& z3)
{
  if (((z2/z1).imag()<=0) || ((z3/z2).imag()<=0)) {
    return false;
  }

  Complex_number<float> one (1,0);
  Complex_number<float> u = (one - conj(z1*z2))   *   (one - conj(z2*z3));
  float a = -(conj(u*z1)*z3).imag();
  float b = (u*(conj(z3-z1))).imag();
  float c = u.imag();

  const float COMPUTATION_TRESHOLD = 0.00001f;
  if (a+b+c> 0 - COMPUTATION_TRESHOLD) {
    return false;
  }

  z0.real(2*c/(std::sqrt(b*b-4*a*c)-b));
  z0.imag(0);
  return true;
}

template<class Traits>
bool
Hyperbolic_fundamental_domain_factory_2<Traits>::
try_to_compute_exact_z3_from_z0_z1_z2(_Cmplx& z0, _Cmplx& z1, _Cmplx& z2, _Cmplx& z3)
{
  _FT zero_number (0);
  _FT one_number (1);
  if ((z0.real()<=zero_number) || (z1.imag()<=zero_number) || (z2.imag()<=zero_number) || (z3.imag()<=zero_number)) {
    return false;
  }

  if ((norm(z0)>=one_number) || (norm(z1)>=one_number) || (norm(z2)>=one_number) || (norm(z3)>=one_number)) {
    return false;
  }

  if (((z1/z0).imag()<=zero_number) || ((z2/z1).imag()<=zero_number) || ((z3/z2).imag()<=zero_number)) {
    return false;
  }

  _Cmplx one_cmplx (_FT(1), _FT(0));
  _Cmplx two_cmplx(_FT(2), _FT(0));

  _Cmplx f_of_z0 = two_cmplx * z0 / (z0*z0 + one_cmplx);
  _Cmplx f_of_z1 = (z0 + z1) / (z0*z1 + one_cmplx);
  _Cmplx f_of_z2 = (z0 + z2) / (z0*z2 + one_cmplx);
  _Cmplx f_of_z3 = (z0 + z3) / (z0*z3 + one_cmplx);

  _Cmplx intermediate = (one_cmplx - f_of_z0*conj(f_of_z1)) * (one_cmplx - f_of_z1*conj(f_of_z2));
  _FT P_of_zero = intermediate.imag();
  _FT P_of_one = (intermediate * (one_cmplx-f_of_z2*conj(f_of_z3))).imag();

  if (P_of_one == P_of_zero) {
    return false;
  }

  _FT lbda = P_of_zero / (P_of_zero - P_of_one);
  _Cmplx V (lbda*(f_of_z3.real()), lbda*(f_of_z3.imag()));

  if ((V.imag()<=zero_number) || (norm(V)>=one_number) || ((V/f_of_z2).imag()<=zero_number)) {
    return false;
  }

  z3 = (V - z0) / (one_cmplx - z0*V);

  return true;
}

////////////////////////////////////////////////////////////////////////////////

template<class Traits>
bool
Hyperbolic_fundamental_domain_factory_2<Traits>::
sanity_check(_Cmplx& z0, _Cmplx& z1, _Cmplx& z2, _Cmplx& z3)
{
  _FT zero_number(0);
  _FT one_number(1);

  // 1. Check the positions
  if ((z0.imag()!=zero_number) || (z0.real()<=zero_number) || (z1.imag()<=zero_number) || (z2.imag()<=zero_number) || (z3.imag()<=zero_number)) {
    return false;
  }

  if ((norm(z0)>=one_number) || (norm(z1)>=one_number) || (norm(z2)>=one_number) || (norm(z3)>=one_number)) {
    return false;
  }

  if (((z2/z1).imag()<=zero_number) || ((z3/z2).imag()<=zero_number)) {
    return false;
  }

  // 2. Check the area
  _Cmplx one_cmplx (one_number, zero_number);
  _Cmplx Z = (one_cmplx-z0*conj(z1)) * (one_cmplx-z1*conj(z2)) *(one_cmplx-z2*conj(z3)) *(one_cmplx+z3*conj(z0));
  if (Z.imag()!=zero_number) {
    return false;
  }

  return true;
}

} // namespace CGAL

#endif // CGAL_HYPERBOLIC_FUNDAMENTAL_DOMAIN_FACTORY_2_H
