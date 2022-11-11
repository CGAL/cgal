// Copyright (c) 1997-2021
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).
// GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sebastien Loriot
//

#ifndef CGAL_INTERSECTIONS_3_BBOX_3_POINT_3_H
#define CGAL_INTERSECTIONS_3_BBOX_3_POINT_3_H

#include <CGAL/Intersection_traits_3.h>
#include <CGAL/Intersections_3/Iso_cuboid_3_Point_3.h>

#include <CGAL/Bbox_3.h>
#include <CGAL/Iso_cuboid_3.h>
#include <CGAL/Point_3.h>

namespace CGAL {

template<typename K>
bool do_intersect(const CGAL::Bbox_3& box,
                  const Point_3<K>& p)
{
  Point_3<K> bl(box.xmin(), box.ymin(), box.zmin()),
             tr(box.xmax(), box.ymax(), box.zmax());
  Iso_cuboid_3<K> ic(bl, tr);
  return K().do_intersect_3_object()(ic, p);
}

template<typename K>
bool do_intersect(const Point_3<K>& a,
                  const CGAL::Bbox_3& b)
{
  return do_intersect(b,a);
}

template<typename K>
typename Intersection_traits<K, typename K::Point_3, CGAL::Bbox_3>::result_type
intersection(const Point_3<K>& p,
             const CGAL::Bbox_3& box)
{
  if(do_intersect(p, box))
    return Intersections::internal::intersection_return<typename K::Intersect_3, typename K::Point_3, CGAL::Bbox_3>(p);

  return Intersections::internal::intersection_return<typename K::Intersect_3, typename K::Point_3, CGAL::Bbox_3>();
}

template<typename K>
typename Intersection_traits<K, CGAL::Bbox_3, typename K::Point_3>::result_type
intersection(const CGAL::Bbox_3& box,
             const Point_3<K>& p)
{
  if(do_intersect(box, p))
    return Intersections::internal::intersection_return<typename K::Intersect_3, CGAL::Bbox_3, typename K::Point_3>(p);

  return Intersections::internal::intersection_return<typename K::Intersect_3, CGAL::Bbox_3, typename K::Point_3>();
}

} // namespace CGAL

#endif // CGAL_INTERSECTIONS_3_BBOX_3_POINT_3_H
