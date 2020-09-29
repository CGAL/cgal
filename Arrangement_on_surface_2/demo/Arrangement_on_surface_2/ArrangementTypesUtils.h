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

#ifndef ARRANGEMENT_DEMO_TYPE_UTILS
#define ARRANGEMENT_DEMO_TYPE_UTILS

#include "ArrangementTypes.h"
#include <CGAL/assertions.h>

namespace demo_types
{

enum class TraitsType : int
{
  SEGMENT_TRAITS,
  POLYLINE_TRAITS,
  LINEAR_TRAITS,
#ifdef CGAL_USE_CORE
  CONIC_TRAITS,
  ALGEBRAIC_TRAITS,
  BEZIER_TRAITS,
  RATIONAL_FUNCTION_TRAITS,
#endif
  NONE,
};

template <class T>
struct TypeHolder
{
  using type = T;
};

namespace details
{
template <typename T>
struct EnumFromTraits
{
  static constexpr TraitsType value = TraitsType::NONE;
};
template <typename Kernel_>
struct EnumFromTraits<CGAL::Arr_segment_traits_2<Kernel_>>
{
  static constexpr TraitsType value = TraitsType::SEGMENT_TRAITS;
};
template <typename SegmentTraits_2_>
struct EnumFromTraits<CGAL::Arr_polyline_traits_2<SegmentTraits_2_>>
{
  static constexpr TraitsType value = TraitsType::POLYLINE_TRAITS;
};
template <typename Kernel_>
struct EnumFromTraits<CGAL::Arr_linear_traits_2<Kernel_>>
{
  static constexpr TraitsType value = TraitsType::LINEAR_TRAITS;
};
#ifdef CGAL_USE_CORE
template <class Rat_kernel_, class Alg_kernel_, class Nt_traits_>
struct EnumFromTraits<
  CGAL::Arr_conic_traits_2<Rat_kernel_, Alg_kernel_, Nt_traits_>>
{
  static constexpr TraitsType value = TraitsType::CONIC_TRAITS;
};
template <
  typename RatKernel_, typename AlgKernel_, typename NtTraits_,
  typename BoundingTraits_>
struct EnumFromTraits<CGAL::Arr_Bezier_curve_traits_2<
  RatKernel_, AlgKernel_, NtTraits_, BoundingTraits_>>
{
  static constexpr TraitsType value = TraitsType::BEZIER_TRAITS;
};
template <typename Coefficient_>
struct EnumFromTraits<CGAL::Arr_algebraic_segment_traits_2<Coefficient_>>
{
  static constexpr TraitsType value = TraitsType::ALGEBRAIC_TRAITS;
};
template <typename AlgebraicKernel_d_1_>
struct EnumFromTraits<
  CGAL::Arr_rational_function_traits_2<AlgebraicKernel_d_1_>>
{
  static constexpr TraitsType value = TraitsType::RATIONAL_FUNCTION_TRAITS;
};
#endif // CGAL_USE_CORE
}

template <typename T>
static constexpr TraitsType enumFromArrType()
{
  return details::EnumFromTraits<typename T::Geometry_traits_2>::value;
}

template <class Lambda>
static void visitArrangementType(TraitsType tt, Lambda lambda)
{
  switch (tt)
  {
  case TraitsType::SEGMENT_TRAITS:
    lambda(TypeHolder<Seg_arr>{});
    break;
  case TraitsType::POLYLINE_TRAITS:
    lambda(TypeHolder<Pol_arr>{});
    break;
  case TraitsType::LINEAR_TRAITS:
    lambda(TypeHolder<Lin_arr>{});
    break;
#ifdef CGAL_USE_CORE
  case TraitsType::CONIC_TRAITS:
    lambda(TypeHolder<Conic_arr>{});
    break;
  case TraitsType::ALGEBRAIC_TRAITS:
    lambda(TypeHolder<Alg_seg_arr>{});
    break;
  case TraitsType::BEZIER_TRAITS:
    lambda(TypeHolder<Bezier_arr>{});
    break;
  case TraitsType::RATIONAL_FUNCTION_TRAITS:
    lambda(TypeHolder<Rational_arr>{});
    break;
#endif
  default:
    CGAL_error();
  }
}

template <class Lambda>
static void forEachArrangementType(Lambda lambda)
{
  lambda(TypeHolder<Seg_arr>{});
  lambda(TypeHolder<Pol_arr>{});
  lambda(TypeHolder<Lin_arr>{});
#ifdef CGAL_USE_CORE
  lambda(TypeHolder<Conic_arr>{});
  lambda(TypeHolder<Alg_seg_arr>{});
  lambda(TypeHolder<Bezier_arr>{});
  lambda(TypeHolder<Rational_arr>{});
#endif
}
}

#endif
