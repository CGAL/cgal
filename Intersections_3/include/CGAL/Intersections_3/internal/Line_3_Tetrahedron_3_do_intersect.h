// Copyright (c) 2019 GeometryFactory(France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri,
//                 Laurent Rineau

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_LINE_3_TETRAHEDRON_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_3_LINE_3_TETRAHEDRON_3_DO_INTERSECT_H

#include <CGAL/Intersections_3/internal/Tetrahedron_3_Unbounded_3_do_intersect.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template<typename K>
typename K::Boolean
do_intersect(const typename K::Line_3& unb,
             const typename K::Tetrahedron_3& tet,
             const K& k)
{
  return do_intersect_tetrahedron_unbounded(tet, unb, k);
}

template<typename K>
typename K::Boolean
do_intersect(const typename K::Tetrahedron_3& tet,
             const typename K::Line_3& unb,
             const K& k)
{
  return do_intersect_tetrahedron_unbounded(tet, unb, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_LINE_3_TETRAHEDRON_3_DO_INTERSECT_H
