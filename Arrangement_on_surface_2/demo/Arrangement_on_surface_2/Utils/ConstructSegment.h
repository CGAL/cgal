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

#ifndef ARRANGEMENT_DEMO_CONSTRUCT_SEGMENT
#define ARRANGEMENT_DEMO_CONSTRUCT_SEGMENT

#include "CurveInputMethods.h"
#include "ForwardDeclarations.h"

// creates a line segment between two points of the type resulting from
// user input (CurveInputMethod::Point_2)
template <typename Traits_>
class Construct_segment
{
public:
  using Traits = Traits_;
  using X_monotone_curve_2 = typename Traits::X_monotone_curve_2;
  using Point_2 = CGAL::Qt::CurveInputMethod::Point_2;

  X_monotone_curve_2
  operator()(const Traits* traits, const Point_2& p1, const Point_2& p2);
};

template <typename Coefficient_>
class Construct_segment<CGAL::Arr_algebraic_segment_traits_2<Coefficient_>>
{
public:
  using Traits = CGAL::Arr_algebraic_segment_traits_2<Coefficient_>;
  using X_monotone_curve_2 = typename Traits::X_monotone_curve_2;
  using Point_2 = CGAL::Qt::CurveInputMethod::Point_2;

  X_monotone_curve_2
  operator()(const Traits* traits, const Point_2& p1, const Point_2& p2);
};

template <
  typename RatKernel, typename AlgKernel, typename NtTraits,
  typename BoundingTraits>
class Construct_segment<CGAL::Arr_Bezier_curve_traits_2<
  RatKernel, AlgKernel, NtTraits, BoundingTraits>>
{
public:
  using Traits = CGAL::Arr_Bezier_curve_traits_2<
    RatKernel, AlgKernel, NtTraits, BoundingTraits>;
  using Curve_2 = typename Traits::Curve_2;
  using X_monotone_curve_2 = typename Traits::X_monotone_curve_2;
  using Point_2 = CGAL::Qt::CurveInputMethod::Point_2;

  X_monotone_curve_2
  operator()(const Traits* traits, const Point_2& p1, const Point_2& p2);
};

template <typename AlgebraicKernel_d_1>
class Construct_segment<
  CGAL::Arr_rational_function_traits_2<AlgebraicKernel_d_1>>
{
public:
  using Traits = CGAL::Arr_rational_function_traits_2<AlgebraicKernel_d_1>;
  using Rational = typename Traits::Rational;
  using Integer = typename Traits::Integer;
  using Polynomial_1 = typename Traits::Polynomial_1;
  using Algebraic_real_1 = typename Traits::Algebraic_real_1;
  using RationalTraits = CGAL::Rational_traits<typename Traits::Rational>;
  using X_monotone_curve_2 = typename Traits::X_monotone_curve_2;
  using Point_2 = CGAL::Qt::CurveInputMethod::Point_2;

  X_monotone_curve_2
  operator()(const Traits* traits, const Point_2& p1, const Point_2& p2);
};

#endif
