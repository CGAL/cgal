// Copyright (c) 2018 INRIA Sophia-Antipolis (France)
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

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_LINE_3_SPHERE_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_3_LINE_3_SPHERE_3_DO_INTERSECT_H

#include <CGAL/Rational_traits.h>
#include <CGAL/Distance_3/Point_3_Line_3.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
inline
typename K::Boolean
do_intersect_impl(const typename K::Line_3& lin,
                  const typename K::Point_3& c,
                  const typename K::FT& r2,
                  const K& k)
{
  typedef typename K::RT RT;
  RT num, den;

  CGAL::internal::squared_distance_RT(c, lin, num, den, k);
  return !(compare_quotients<RT>(num, den,
                                 Rational_traits<typename K::FT>().numerator(r2),
                                 Rational_traits<typename K::FT>().denominator(r2)) == LARGER);
}

template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Sphere_3& sp,
             const typename K::Line_3& lin,
             const K& k)
{
  return do_intersect_impl(lin, sp.center(), sp.squared_radius(), k);
}

template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Line_3& lin,
             const typename K::Sphere_3& sp,
             const K& k)
{
  return do_intersect(sp, lin, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_LINE_3_SPHERE_3_DO_INTERSECT_H
