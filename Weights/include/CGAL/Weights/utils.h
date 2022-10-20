// Copyright (c) 2020 GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Dmitry Anisimov
//

#ifndef CGAL_WEIGHTS_UTILS_H
#define CGAL_WEIGHTS_UTILS_H

#include <CGAL/Weights/internal/utils.h>

#include <boost/algorithm/clamp.hpp>

namespace CGAL {
namespace Weights {

/// \cond SKIP_IN_MANUAL

// Computes cotangent between two 2D vectors.
template<typename GeomTraits>
typename GeomTraits::FT cotangent_2(const typename GeomTraits::Point_2& p,
                                    const typename GeomTraits::Point_2& q,
                                    const typename GeomTraits::Point_2& r,
                                    const GeomTraits& traits)
{
  using FT = typename GeomTraits::FT;
  using Vector_2 = typename GeomTraits::Vector_2;

  auto dot_product_2 = traits.compute_scalar_product_2_object();
  auto determinant_2 = traits.compute_determinant_2_object();
  auto vector_2 = traits.construct_vector_2_object();

  const Vector_2 v1 = vector_2(q, r);
  const Vector_2 v2 = vector_2(q, p);

  const FT dot = dot_product_2(v1, v2);
  const FT length = CGAL::abs(determinant_2(v1, v2));

  if (!is_zero(length))
    return dot / length;

  return FT(0); // undefined
}

template<typename GeomTraits>
typename GeomTraits::FT cotangent(const typename GeomTraits::Point_2& p,
                                  const typename GeomTraits::Point_2& q,
                                  const typename GeomTraits::Point_2& r,
                                  const GeomTraits& traits)
{
  return cotangent_2(p, q, r, traits);
}

template<typename Kernel>
typename Kernel::FT cotangent(const CGAL::Point_2<Kernel>& p,
                              const CGAL::Point_2<Kernel>& q,
                              const CGAL::Point_2<Kernel>& r)
{
  const Kernel traits;
  return cotangent(p, q, r, traits);
}

// =================================================================================================

// Computes cotangent between two 3D vectors.
template<typename GeomTraits>
typename GeomTraits::FT cotangent_3(const typename GeomTraits::Point_3& p,
                                    const typename GeomTraits::Point_3& q,
                                    const typename GeomTraits::Point_3& r,
                                    const GeomTraits& traits)
{
  using FT = typename GeomTraits::FT;
  using Vector_3 = typename GeomTraits::Vector_3;

  auto dot_product_3 = traits.compute_scalar_product_3_object();
  auto cross_product_3 = traits.construct_cross_product_vector_3_object();
  auto vector_3 = traits.construct_vector_3_object();

  const Vector_3 v1 = vector_3(q, r);
  const Vector_3 v2 = vector_3(q, p);

  const FT dot = dot_product_3(v1, v2);
  const Vector_3 cross = cross_product_3(v1, v2);

  const FT length = internal::length_3(cross, traits);
  if (!is_zero(length))
    return dot / length;

  return FT(0); // undefined
}

template<typename GeomTraits>
typename GeomTraits::FT cotangent(const typename GeomTraits::Point_3& p,
                                  const typename GeomTraits::Point_3& q,
                                  const typename GeomTraits::Point_3& r,
                                  const GeomTraits& traits)
{
  return cotangent_3(p, q, r, traits);
}

template<typename Kernel>
typename Kernel::FT cotangent(const CGAL::Point_3<Kernel>& p,
                              const CGAL::Point_3<Kernel>& q,
                              const CGAL::Point_3<Kernel>& r)
{
  const Kernel traits;
  return cotangent(p, q, r, traits);
}

// =================================================================================================

// Computes tangent between two 2D vectors.
template<typename GeomTraits>
typename GeomTraits::FT tangent_2(const typename GeomTraits::Point_2& p,
                                  const typename GeomTraits::Point_2& q,
                                  const typename GeomTraits::Point_2& r,
                                  const GeomTraits& traits)
{
  using FT = typename GeomTraits::FT;
  using Vector_2 = typename GeomTraits::Vector_2;

  auto dot_product_2 = traits.compute_scalar_product_2_object();
  auto determinant_2 = traits.compute_determinant_2_object();
  auto vector_2 = traits.construct_vector_2_object();

  const Vector_2 v1 = vector_2(q, r);
  const Vector_2 v2 = vector_2(q, p);

  const FT dot = dot_product_2(v1, v2);
  if (!is_zero(dot))
  {
    const FT cross = determinant_2(v1, v2);
    const FT length = CGAL::abs(cross);
    return length / dot;
  }

  return FT(0); // undefined
}

template<typename GeomTraits>
typename GeomTraits::FT tangent(const typename GeomTraits::Point_2& p,
                                const typename GeomTraits::Point_2& q,
                                const typename GeomTraits::Point_2& r,
                                const GeomTraits& traits)
{
  return tangent_2(p, q, r, traits);
}

template<typename Kernel>
typename Kernel::FT tangent(const CGAL::Point_2<Kernel>& p,
                            const CGAL::Point_2<Kernel>& q,
                            const CGAL::Point_2<Kernel>& r)
{
  const Kernel traits;
  return tangent(p, q, r, traits);
}

// =================================================================================================

// Computes tangent between two 3D vectors.
template<typename GeomTraits>
typename GeomTraits::FT tangent_3(const typename GeomTraits::Point_3& p,
                                  const typename GeomTraits::Point_3& q,
                                  const typename GeomTraits::Point_3& r,
                                  const GeomTraits& traits)
{
  using FT = typename GeomTraits::FT;
  using Vector_3 = typename GeomTraits::Vector_3;

  auto dot_product_3 = traits.compute_scalar_product_3_object();
  auto cross_product_3 = traits.construct_cross_product_vector_3_object();
  auto vector_3 = traits.construct_vector_3_object();

  const Vector_3 v1 = vector_3(q, r);
  const Vector_3 v2 = vector_3(q, p);

  const FT dot = dot_product_3(v1, v2);
  if (!is_zero(dot))
  {
    const Vector_3 cross = cross_product_3(v1, v2);
    const FT length = internal::length_3(cross, traits);
    return length / dot;
  }

  return FT(0); // undefined
}

template<typename GeomTraits>
typename GeomTraits::FT tangent(const typename GeomTraits::Point_3& p,
                                const typename GeomTraits::Point_3& q,
                                const typename GeomTraits::Point_3& r,
                                const GeomTraits& traits)
{
  return tangent_3(p, q, r, traits);
}

template<typename Kernel>
typename Kernel::FT tangent(const CGAL::Point_3<Kernel>& p,
                            const CGAL::Point_3<Kernel>& q,
                            const CGAL::Point_3<Kernel>& r)
{
  const Kernel traits;
  return tangent(p, q, r, traits);
}

// =================================================================================================

// Computes a clamped cotangent between two 3D vectors.
// In the old version of weights in PMP, it was called "Cotangent_value_Meyer_secure".
// See Weights/internal/pmp_weights_deprecated.h for more information.
template<typename GeomTraits>
typename GeomTraits::FT cotangent_3_clamped(const typename GeomTraits::Point_3& p,
                                            const typename GeomTraits::Point_3& q,
                                            const typename GeomTraits::Point_3& r,
                                            const GeomTraits& traits)
{
  using FT = typename GeomTraits::FT;
  using Vector_3 = typename GeomTraits::Vector_3;

  using Get_sqrt = internal::Get_sqrt<GeomTraits>;
  auto sqrt = Get_sqrt::sqrt_object(traits);

  auto dot_product_3 = traits.compute_scalar_product_3_object();
  auto vector_3 = traits.construct_vector_3_object();

  const Vector_3 v1 = vector_3(q, r);
  const Vector_3 v2 = vector_3(q, p);

  const FT dot = dot_product_3(v1, v2);
  const FT length_v1 = internal::length_3(v1, traits);
  const FT length_v2 = internal::length_3(v2, traits);

  const FT lb = -FT(999) / FT(1000),
           ub =  FT(999) / FT(1000);
  const FT cosine = boost::algorithm::clamp<FT>(dot / (length_v1 * length_v2), lb, ub);
  const FT sine = sqrt(FT(1) - square(cosine));

  CGAL_assertion(!is_zero(sine));
  if (!is_zero(sine))
    return cosine / sine;

  return FT(0); // undefined
}

/// \endcond

} // namespace Weights
} // namespace CGAL

#endif // CGAL_WEIGHTS_UTILS_H
