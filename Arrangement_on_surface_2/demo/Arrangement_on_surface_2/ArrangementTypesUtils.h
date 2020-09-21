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
#include <type_traits>

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

template <typename T>
static constexpr TraitsType enumFromArrType()
{
  using namespace std;

  if (is_same<T, Seg_arr>::value) return TraitsType::SEGMENT_TRAITS;
  else if (is_same<T, Pol_arr>::value) return TraitsType::POLYLINE_TRAITS;
  else if (is_same<T, Lin_arr>::value) return TraitsType::LINEAR_TRAITS;
#ifdef CGAL_USE_CORE
  else if (is_same<T, Conic_arr>::value) return TraitsType::CONIC_TRAITS;
  else if (is_same<T, Alg_seg_arr>::value) return TraitsType::ALGEBRAIC_TRAITS;
  else if (is_same<T, Bezier_arr>::value) return TraitsType::BEZIER_TRAITS;
  else if (is_same<T, Rational_arr>::value) return TraitsType::RATIONAL_FUNCTION_TRAITS;
#endif
  else return TraitsType::NONE;
}

template <class Lambda>
static void visitArrangementType(TraitsType tt, Lambda lambda)
{
  switch (tt)
  {
  default:
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
