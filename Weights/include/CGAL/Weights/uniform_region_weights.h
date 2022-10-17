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

/// \cond SKIP_IN_MANUAL
template<typename GeomTraits>
typename GeomTraits::FT uniform_area(const typename GeomTraits::Point_2&,
                                     const typename GeomTraits::Point_2&,
                                     const typename GeomTraits::Point_2&,
                                     const GeomTraits&)
{
  using FT = typename GeomTraits::FT;
  return FT(1);
}

template<typename GeomTraits>
typename GeomTraits::FT uniform_area(const CGAL::Point_2<GeomTraits>& p,
                                     const CGAL::Point_2<GeomTraits>& q,
                                     const CGAL::Point_2<GeomTraits>& r)
{
  const GeomTraits traits;
  return uniform_area(p, q, r, traits);
}

template<typename GeomTraits>
typename GeomTraits::FT uniform_area(const typename GeomTraits::Point_3&,
                                     const typename GeomTraits::Point_3&,
                                     const typename GeomTraits::Point_3&,
                                     const GeomTraits&)
{
  using FT = typename GeomTraits::FT;
  return FT(1);
}

template<typename GeomTraits>
typename GeomTraits::FT uniform_area(const CGAL::Point_3<GeomTraits>& p,
                                     const CGAL::Point_3<GeomTraits>& q,
                                     const CGAL::Point_3<GeomTraits>& r)
{
  const GeomTraits traits;
  return uniform_area(p, q, r, traits);
}

/// \endcond

} // namespace Weights
} // namespace CGAL

#endif // CGAL_UNIFORM_REGION_WEIGHTS_H
