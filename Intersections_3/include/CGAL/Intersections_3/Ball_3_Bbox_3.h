// Copyright (c) 2026 GeometryFactory Sarl
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

#ifndef CGAL_INTERSECTIONS_3_BALL_3_BBOX_3_H
#define CGAL_INTERSECTIONS_3_BALL_3_BBOX_3_H

#include <CGAL/Intersection_traits_3.h>
#include <CGAL/Intersections_3/internal/Ball_3_Bbox_3_do_intersect.h>

#include <CGAL/Ball_3.h>
#include <CGAL/Bbox_3.h>

namespace CGAL {

template<typename K>
typename K::Boolean
do_intersect(const CGAL::Bbox_3& box,
             const Ball_3<K>& b)
{
  return K().do_intersect_3_object()(box, b);
}

template<typename K>
typename K::Boolean
do_intersect(const Ball_3<K>& b,
             const CGAL::Bbox_3& box)
{
  return K().do_intersect_3_object()(b, box);
}

} // namespace CGAL

#endif // CGAL_INTERSECTIONS_3_BALL_3_BBOX_3_H
