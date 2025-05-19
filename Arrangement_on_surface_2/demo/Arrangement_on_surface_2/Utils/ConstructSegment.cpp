// Copyright (c) 2020 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Ahmed Essam <theartful.ae@gmail.com>

#include "ConstructSegment.h"
#include "ArrangementTypes.h"
#include "Utils.h"

#include <CGAL/polynomial_utils.h>

template <typename Traits_>
auto Construct_segment<Traits_>::operator()(
  const Traits* traits, const Point_2& p1, const Point_2& p2)
  -> X_monotone_curve_2
{
  Arr_construct_point_2<Traits> toArrPoint{traits};

  auto construct_x_monotone_curve_2 =
    traits->construct_x_monotone_curve_2_object();
  return construct_x_monotone_curve_2(toArrPoint(p1), toArrPoint(p2));
}

template <typename Coefficient_>
auto Construct_segment<CGAL::Arr_algebraic_segment_traits_2<Coefficient_>>::
operator()(const Traits* traits, const Point_2& p1, const Point_2& p2)
  -> X_monotone_curve_2
{
  Arr_construct_point_2<Traits> toArrPoint{traits};

  auto construct_x_monotone_segment_2 =
    traits->construct_x_monotone_segment_2_object();

  std::vector<X_monotone_curve_2> curves;
  construct_x_monotone_segment_2(
    toArrPoint(p1), toArrPoint(p2), std::back_inserter(curves));
  return curves[0];
}

template <
  typename RatKernel, typename AlgKernel, typename NtTraits,
  typename BoundingTraits>
auto Construct_segment<CGAL::Arr_Bezier_curve_traits_2<
  RatKernel, AlgKernel, NtTraits, BoundingTraits>>::
operator()(const Traits* traits, const Point_2& p1, const Point_2& p2)
  -> X_monotone_curve_2
{
  Arr_construct_point_2<Traits> toArrPoint{traits};

  typename Traits::Point_2 points[] = {toArrPoint(p1), toArrPoint(p2)};
  Curve_2 curve{points, points + 2};
  auto make_x_monotone = traits->make_x_monotone_2_object();
  std::vector<CGAL::Object> curves;
  make_x_monotone(curve, std::back_inserter(curves));
  X_monotone_curve_2 segment;
  CGAL::assign(segment, curves[0]);
  return segment;
}

template <typename AlgebraicKernel_d_1>
auto Construct_segment<
  CGAL::Arr_rational_function_traits_2<AlgebraicKernel_d_1>>::
operator()(const Traits* traits, const Point_2& p1, const Point_2& p2)
  -> X_monotone_curve_2
{
  RationalTraits ratTraits;

  Rational dx = p2.x() - p1.x();
  Rational dy = p2.y() - p1.y();
  Polynomial_1 x = CGAL::shift(Polynomial_1(1), 1, 0);
  Polynomial_1 poly_num;
  Polynomial_1 poly_den;

  if (dx != 0)
  {
    Rational mRat = dy / dx;
    Rational cRat = p1.y() - mRat * p1.x();
    // y = (a/b) x + (e/f)
    auto a = ratTraits.numerator(mRat);
    auto b = ratTraits.denominator(mRat);
    auto e = ratTraits.numerator(cRat);
    auto f = ratTraits.denominator(cRat);

    poly_num = f * a * x + b * e;
    poly_den = Polynomial_1{b * f};
  }
  else
  {
    throw std::runtime_error(
      "Vertical split segments are not allowed with rational traits!\n");
  }

  auto construct_x_monotone_curve_2 =
    traits->construct_x_monotone_curve_2_object();

  return construct_x_monotone_curve_2(
    poly_num, poly_den, Algebraic_real_1{p1.x()}, Algebraic_real_1{p2.x()});
}

ARRANGEMENT_DEMO_SPECIALIZE_TRAITS(Construct_segment)
