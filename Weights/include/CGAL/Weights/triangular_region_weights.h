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

/// \cond SKIP_IN_MANUAL
template<typename GeomTraits>
typename GeomTraits::FT triangular_area(const typename GeomTraits::Point_2& p,
                                        const typename GeomTraits::Point_2& q,
                                        const typename GeomTraits::Point_2& r,
                                        const GeomTraits& traits)
{
  return internal::positive_area_2(traits, p, q, r);
}

template<typename GeomTraits>
typename GeomTraits::FT triangular_area(const CGAL::Point_2<GeomTraits>& p,
                                        const CGAL::Point_2<GeomTraits>& q,
                                        const CGAL::Point_2<GeomTraits>& r)
{
  const GeomTraits traits;
  return triangular_area(p, q, r, traits);
}

template<typename GeomTraits>
typename GeomTraits::FT triangular_area(const typename GeomTraits::Point_3& p,
                                        const typename GeomTraits::Point_3& q,
                                        const typename GeomTraits::Point_3& r,
                                        const GeomTraits& traits)
{
  return internal::positive_area_3(traits, p, q, r);
}

template<typename GeomTraits>
typename GeomTraits::FT triangular_area(const CGAL::Point_3<GeomTraits>& p,
                                        const CGAL::Point_3<GeomTraits>& q,
                                        const CGAL::Point_3<GeomTraits>& r)
{
  const GeomTraits traits;
  return triangular_area(p, q, r, traits);
}

/// \endcond

} // namespace Weights
} // namespace CGAL

#endif // CGAL_TRIANGULAR_REGION_WEIGHTS_H
