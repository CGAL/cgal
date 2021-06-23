// Copyright (c) 2018 GeometryFactory Sarl
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_INTERNAL_INTERSECTIONS_SPHERE_3_TRIANGLE_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_SPHERE_3_TRIANGLE_3_DO_INTERSECT_H

#include <CGAL/Rational_traits.h>
#include <CGAL/squared_distance_3_2.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Sphere_3& sp,
             const typename K::Triangle_3& tr,
             const K& k)
{
  typedef typename K::RT RT;
  RT num, den;

  CGAL::internal::squared_distance_RT(sp.center(), tr, num, den, k);
  return !(compare_quotients<RT>(num, den,
                                 Rational_traits<typename K::FT>().numerator(sp.squared_radius()),
                                 Rational_traits<typename K::FT>().denominator(sp.squared_radius())) == LARGER);
}

template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Triangle_3& tr,
             const typename K::Sphere_3& sp,
             const K& k)
{
  return do_intersect(sp, tr, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_SPHERE_3_TRIANGLE_3_DO_INTERSECT_H
