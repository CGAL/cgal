// Copyright (c) 2008  INRIA Sophia-Antipolis (France), ETH Zurich (Switzerland).
// Copyright (c) 2010, 2014  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Camille Wormser, Jane Tournois, Pierre Alliez

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_ISO_CUBOID_3_TRIANGLE_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_3_ISO_CUBOID_3_TRIANGLE_3_DO_INTERSECT_H

#include <CGAL/Intersections_3/internal/Bbox_3_Triangle_3_do_intersect.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
bool do_intersect(const typename K::Triangle_3& triangle,
                  const typename K::Iso_cuboid_3& bbox,
                  const K& k)
{
  return do_intersect_bbox_or_iso_cuboid(triangle, bbox, k);
}

template <class K>
bool do_intersect(const typename K::Iso_cuboid_3& bbox,
                  const typename K::Triangle_3& triangle,
                  const K& k)
{
  return do_intersect_bbox_or_iso_cuboid(triangle, bbox, k);
}

} // namespace internal
} // namespace Intersections



} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_ISO_CUBOID_3_TRIANGLE_3_DO_INTERSECT_H
