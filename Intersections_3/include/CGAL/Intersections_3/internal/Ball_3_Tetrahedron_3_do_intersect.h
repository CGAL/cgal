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

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_BALL_3_TETRAHEDRON_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_3_BALL_3_TETRAHEDRON_3_DO_INTERSECT_H

#include <CGAL/Intersections_3/internal/Ball_3_Triangle_3_do_intersect.h>
#include <CGAL/Intersections_3/internal/Tetrahedron_3_Bounded_3_do_intersect.h>

namespace CGAL {
namespace Intersections {
namespace internal {

// The ball intersects one of the faces of the tetrahedron (as a solid), or
// it is strictly inside the tetrahedron (then its center is inside).
template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Ball_3& b,
             const typename K::Tetrahedron_3& tet,
             const K& k)
{
  return do_intersect_tetrahedron_bounded(b, tet, b.center(), k);
}

template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Tetrahedron_3& tet,
             const typename K::Ball_3& b,
             const K& k)
{
  return do_intersect(b, tet, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_BALL_3_TETRAHEDRON_3_DO_INTERSECT_H
