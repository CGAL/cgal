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


#ifndef CGAL_INTERSECTIONS_BBOX_2_ISO_RECTANGLE_2_H
#define CGAL_INTERSECTIONS_BBOX_2_ISO_RECTANGLE_2_H

#include <CGAL/Bbox_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Intersections_2/Iso_rectangle_2_Iso_rectangle_2.h>

namespace CGAL {

template <typename K>
inline bool do_intersect(const Iso_rectangle_2<K> &rect,
                         const Bbox_2 &box)
{
  return do_intersect(K::Iso_rectangle_2(box), rect);
}

template <typename K>
inline bool do_intersect(const Bbox_2 &box,
                         const Iso_rectangle_2<K> &rect)
{
  return do_intersect(rect, box);
}

template<typename K>
typename Intersection_traits<K, typename K::Iso_rectangle_2, Bbox_2>::result_type
intersection(const Bbox_2& box,
             const Iso_rectangle_2<K>& rect)
{
  typename K::Iso_rectangle_2 rec(box.xmin(), box.ymin(), box.xmax(), box.ymax());
  return intersection(rec, rect);
}

template<typename K>
typename Intersection_traits<K, typename K::Iso_rectangle_2, Bbox_2>::result_type
intersection(const Iso_rectangle_2<K>& rect,
             const Bbox_2& box)
{
  return intersection(box, rect);
}

} // namespace CGAL

#endif // CGAL_INTERSECTIONS_BBOX_2_ISO_RECTANGLE_2_H
