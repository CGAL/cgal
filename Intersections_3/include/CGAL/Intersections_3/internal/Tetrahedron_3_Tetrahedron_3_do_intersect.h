// Copyright (c) 2005
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Nico Kruithof

#ifndef CGAL_INTERNAL_INTERSECTIONS_TETRAHEDRON_3_TETRAHEDRON_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_TETRAHEDRON_3_TETRAHEDRON_3_DO_INTERSECT_H

#include <CGAL/Intersections_3/internal/Tetrahedron_3_Triangle_3_do_intersect.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Tetrahedron_3& lh_tet,
             const typename K::Tetrahedron_3& rh_tet,
             const K& k)
{
  return do_intersect_tetrahedron_bounded(lh_tet, rh_tet, lh_tet[0], k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_TETRAHEDRON_3_TETRAHEDRON_3_DO_INTERSECT_H
