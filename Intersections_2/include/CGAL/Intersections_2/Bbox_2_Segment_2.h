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


#ifndef CGAL_INTERSECTIONS_BBOX_2_SEGMENT_2_H
#define CGAL_INTERSECTIONS_BBOX_2_SEGMENT_2_H

#include <CGAL/Bbox_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Intersections_2/Iso_rectangle_2_Segment_2.h>

namespace CGAL {


template <class K>
inline bool do_intersect(
    const Segment_2<K> &seg,
    const Bbox_2 &box)
{
  typename K::Iso_rectangle_2 rec(box.xmin(), box.ymin(), box.xmax(), box.ymax());
  return do_intersect(rec, seg);
}

template <class K>
inline bool do_intersect(
    const Bbox_2 &box,
    const Segment_2<K> &seg)
{
  return do_intersect(seg, box);
}

template<typename K>
typename Intersection_traits<K, typename K::Segment_2, Bbox_2>::result_type
intersection(const CGAL::Bbox_2& box,
             const Segment_2<K>& seg) {
  typename K::Iso_rectangle_2 rec(box.xmin(), box.ymin(), box.xmax(), box.ymax());
  return intersection(rec, seg);
}

template<typename K>
typename Intersection_traits<K, typename K::Segment_2, Bbox_2>::result_type
intersection(const Segment_2<K>& seg,
             const CGAL::Bbox_2& box) {
  return intersection(box, seg);
}

}
#endif // CGAL_INTERSECTIONS_BBOX_2_SEGMENT_2_H
