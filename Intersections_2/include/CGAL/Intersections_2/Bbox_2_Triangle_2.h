// Copyright (c) 2019
// GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
//
// Author(s)     : Maxime Gimeno


#ifndef CGAL_INTERSECTIONS_BBOX_2_TRIANGLE_2_H
#define CGAL_INTERSECTIONS_BBOX_2_TRIANGLE_2_H

#include <CGAL/Bbox_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/kernel_assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/Intersections_2/Iso_rectangle_2_Triangle_2.h>

namespace CGAL {


template <class K>
inline bool do_intersect(
    const Triangle_2<K> &line,
    const Bbox_2 &box)
{
  typename K::Iso_rectangle_2 rec(box.xmin(), box.ymin(), box.xmax(), box.ymax());
  return do_intersect(rec, line);
}

template <class K>
inline bool do_intersect(
    const Bbox_2 &box,
const Triangle_2<K> &line)
{
  return do_intersect(line, box);
}

template<typename K>
typename Intersection_traits<K, typename K::Triangle_2, Bbox_2>::result_type
intersection(const CGAL::Bbox_2& box,
             const Triangle_2<K>& line) {
  typename K::Iso_rectangle_2 rec(box.xmin(), box.ymin(), box.xmax(), box.ymax());
  return intersection(rec, line);
}

template<typename K>
typename Intersection_traits<K, typename K::Triangle_2, Bbox_2>::result_type
intersection(const Triangle_2<K>& line,
             const CGAL::Bbox_2& box) {
  return intersection(box, line);
}

}
#endif // BBOX_2_TRIANGLE_2_H
