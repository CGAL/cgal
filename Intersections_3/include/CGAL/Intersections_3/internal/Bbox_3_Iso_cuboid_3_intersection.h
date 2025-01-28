// Copyright (c) 1997-2010
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
// Author(s)     : Geert-Jan Giezeman <geert@cs.uu.nl>
//                 Sebastien Loriot <Sebastien.Loriot@geometryfactory.com>

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_ISO_CUBOID_3_INTERSECTION_H
#define CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_ISO_CUBOID_3_INTERSECTION_H

#include <CGAL/Intersection_traits_3.h>
#include <CGAL/Intersections_3/internal/Iso_cuboid_3_Iso_cuboid_3_intersection.h>

#include <CGAL/Bbox_3.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
typename Intersection_traits<K, CGAL::Bbox_3, typename K::Iso_cuboid_3>::result_type
intersection(const CGAL::Bbox_3& box,
             const typename K::Iso_cuboid_3& cub,
             const K&)
{
  typename K::Iso_cuboid_3 iso_cub(box);
  return intersection(iso_cub, cub);
}

template <class K>
inline
typename Intersection_traits<K, typename K::Iso_cuboid_3, CGAL::Bbox_3>::result_type
intersection(const typename K::Iso_cuboid_3& cub,
             const CGAL::Bbox_3& box,
             const K& k)
{
  return intersection(box, cub, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_ISO_CUBOID_3_INTERSECTION_H
