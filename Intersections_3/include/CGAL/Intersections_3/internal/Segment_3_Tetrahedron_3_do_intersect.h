// Copyright (c) 2001 GeometryFactory(France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Nico Kruithof

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_SEGMENT_3_TETRAHEDRON_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_3_SEGMENT_3_TETRAHEDRON_3_DO_INTERSECT_H

#include <CGAL/Intersections_3/internal/Tetrahedron_3_Bounded_3_do_intersect.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Tetrahedron_3& tet,
             const typename K::Segment_3& seg,
             const K& k)
{
  return do_intersect_tetrahedron_bounded(seg, tet, seg.source(), k);
}

template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Segment_3& seg,
             const typename K::Tetrahedron_3& tet,
             const K& k)
{
  return do_intersect_tetrahedron_bounded(seg, tet, seg.source(), k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_SEGMENT_3_TETRAHEDRON_3_DO_INTERSECT_H
