// Copyright (c) 2026 GeometryFactory Sarl
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sébastien Loriot
//

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_BALL_3_LINE_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_3_BALL_3_LINE_3_DO_INTERSECT_H

#include <CGAL/Rational_traits.h>
#include <CGAL/Intersections_3/internal/Line_3_Sphere_3_do_intersect.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Line_3& l,
             const typename K::Ball_3& b,
             const K& k)
{
  return do_intersect_impl(l, b.center(), b.squared_radius(), k);
}

template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Ball_3& b,
             const typename K::Line_3& l,
             const K& k)
{
  return do_intersect(l, b, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_BALL_3_LINE_3_DO_INTERSECT_H
