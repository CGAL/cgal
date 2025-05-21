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

#ifndef CGAL_UNIFORM_REGION_WEIGHTS_H
#define CGAL_UNIFORM_REGION_WEIGHTS_H

#include <CGAL/Point_2.h>
#include <CGAL/Point_3.h>

namespace CGAL {
namespace Weights {

// 2D ==============================================================================================

/*!
  \ingroup PkgWeightsRefUniformRegionWeights
  \brief this function always returns `1`, given three 2D points.
  \tparam GeomTraits a model of `AnalyticWeightTraits_2`
*/
template<typename GeomTraits>
typename GeomTraits::FT uniform_area(const typename GeomTraits::Point_2&,
                                     const typename GeomTraits::Point_2&,
                                     const typename GeomTraits::Point_2&,
                                     const GeomTraits&)
{
  using FT = typename GeomTraits::FT;
  return FT(1);
}

/*!
  \ingroup PkgWeightsRefUniformRegionWeights
  \brief this function always returns `1`, given three 2D points in 2D.
  \tparam Kernel a model of `Kernel`
*/
template<typename Kernel>
typename Kernel::FT uniform_area(const CGAL::Point_2<Kernel>& p,
                                 const CGAL::Point_2<Kernel>& q,
                                 const CGAL::Point_2<Kernel>& r)
{
  const Kernel traits;
  return uniform_area(p, q, r, traits);
}

// 3D ==============================================================================================

/*!
  \ingroup PkgWeightsRefUniformRegionWeights
  \brief this function always returns `1`, given three 3D points.
  \tparam GeomTraits a model of `AnalyticWeightTraits_3`
*/
template<typename GeomTraits>
typename GeomTraits::FT uniform_area(const typename GeomTraits::Point_3&,
                                     const typename GeomTraits::Point_3&,
                                     const typename GeomTraits::Point_3&,
                                     const GeomTraits&)
{
  using FT = typename GeomTraits::FT;
  return FT(1);
}

/*!
  \ingroup PkgWeightsRefUniformRegionWeights
  \brief this function always returns `1`, given three 3D points.
  \tparam Kernel a model of `Kernel`
*/
template<typename Kernel>
typename Kernel::FT uniform_area(const CGAL::Point_3<Kernel>& p,
                                 const CGAL::Point_3<Kernel>& q,
                                 const CGAL::Point_3<Kernel>& r)
{
  const Kernel traits;
  return uniform_area(p, q, r, traits);
}

} // namespace Weights
} // namespace CGAL

#endif // CGAL_UNIFORM_REGION_WEIGHTS_H
