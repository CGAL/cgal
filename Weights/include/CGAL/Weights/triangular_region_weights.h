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

#ifndef CGAL_TRIANGULAR_REGION_WEIGHTS_H
#define CGAL_TRIANGULAR_REGION_WEIGHTS_H

#include <CGAL/Weights/internal/utils.h>

#include <CGAL/Point_2.h>
#include <CGAL/Point_3.h>

namespace CGAL {
namespace Weights {

// 2D ==============================================================================================

/*!
  \ingroup PkgWeightsRefTriangularRegionWeights
  \brief computes the area of the triangular cell in 2D using the points `p`, `q`, and `r`
  \tparam GeomTraits a model of `AnalyticWeightTraits_2`
*/
template<typename GeomTraits>
typename GeomTraits::FT triangular_area(const typename GeomTraits::Point_2& p,
                                        const typename GeomTraits::Point_2& q,
                                        const typename GeomTraits::Point_2& r,
                                        const GeomTraits& traits)
{
  return internal::positive_area_2(p, q, r, traits);
}

/*!
  \ingroup PkgWeightsRefTriangularRegionWeights
  \brief computes the area of the triangular cell in 2D using the points `p`, `q`, and `r`.
  \tparam Kernel a model of `Kernel`
*/
template<typename Kernel>
typename Kernel::FT triangular_area(const CGAL::Point_2<Kernel>& p,
                                    const CGAL::Point_2<Kernel>& q,
                                    const CGAL::Point_2<Kernel>& r)
{
  const Kernel traits;
  return triangular_area(p, q, r, traits);
}

// 3D ==============================================================================================

/*!
  \ingroup PkgWeightsRefTriangularRegionWeights
  \brief computes the area of the triangular cell in 3D using the points `p`, `q`, and `r`.
  \tparam GeomTraits a model of `AnalyticWeightTraits_3`
*/
template<typename GeomTraits>
typename GeomTraits::FT triangular_area(const typename GeomTraits::Point_3& p,
                                        const typename GeomTraits::Point_3& q,
                                        const typename GeomTraits::Point_3& r,
                                        const GeomTraits& traits)
{
  return internal::positive_area_3(p, q, r, traits);
}

/*!
  \ingroup PkgWeightsRefTriangularRegionWeights
  \brief computes the area of the triangular cell in 3D using the points `p`, `q`, and `r`.
  \tparam Kernel a model of `Kernel`
*/
template<typename Kernel>
typename Kernel::FT triangular_area(const CGAL::Point_3<Kernel>& p,
                                    const CGAL::Point_3<Kernel>& q,
                                    const CGAL::Point_3<Kernel>& r)
{
  const Kernel traits;
  return triangular_area(p, q, r, traits);
}

} // namespace Weights
} // namespace CGAL

#endif // CGAL_TRIANGULAR_REGION_WEIGHTS_H
