// Copyright (c) 2019
// GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Maxime Gimeno


#ifndef CGAL_INTERSECTIONS_BBOX_2_TRIANGLE_2_H
#define CGAL_INTERSECTIONS_BBOX_2_TRIANGLE_2_H

#include <CGAL/Bbox_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Intersections_2/Iso_rectangle_2_Triangle_2.h>

namespace CGAL {


template <class K>
inline bool do_intersect(
    const Triangle_2<K> &tr,
    const Bbox_2 &box)
{
  typename K::Iso_rectangle_2 rec(box.xmin(), box.ymin(), box.xmax(), box.ymax());
  return do_intersect(rec, tr);
}

template <class K>
inline bool do_intersect(
    const Bbox_2 &box,
    const Triangle_2<K> &tr)
{
  return do_intersect(tr, box);
}

template<typename K>
typename Intersection_traits<K, typename K::Triangle_2, Bbox_2>::result_type
intersection(const Bbox_2& box,
             const Triangle_2<K>& tr) {
  typename K::Iso_rectangle_2 rec(box.xmin(), box.ymin(), box.xmax(), box.ymax());
  return intersection(rec, tr);
}

template<typename K>
typename Intersection_traits<K, typename K::Triangle_2, Bbox_2>::result_type
intersection(const Triangle_2<K>& tr,
             const Bbox_2& box) {
  return intersection(box, tr);
}

}
#endif // CGAL_INTERSECTIONS_BBOX_2_TRIANGLE_2_H
