// Copyright (c) 2019
// GeometryFactory.  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Maxime Gimeno

#ifndef CGAL_INTERSECTIONS_2_CIRCLE_2_TRIANGLE_2_H
#define CGAL_INTERSECTIONS_2_CIRCLE_2_TRIANGLE_2_H

#include <CGAL/Circle_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/squared_distance_2_1.h>
#include <CGAL/Intersection_traits_2.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
bool
do_intersect(const typename K::Circle_2 & c,
             const typename K::Triangle_2& t,
             const K&)
{
  char nb_close(0);
//if all three points are inside the circle, there is no intersection.
  if(squared_distance(c.center(), t) < c.squared_radius())
  {
    for(int i = 0; i< 3; ++i)
    {
      if(squared_distance(c.center(), t[i]) == c.squared_radius())
        return true;
      if(squared_distance(c.center(), t[i]) < c.squared_radius())
        ++nb_close;
    }
    return (nb_close < 3);
  }
  else
    return squared_distance(c.center(), t) == c.squared_radius();
}

template <class K>
bool
do_intersect(const typename K::Triangle_2& t,
             const typename K::Circle_2 & c,
             const K&)
{
  return do_intersect(c,t);
}

} // namespace internal
} // namespace Intersections

CGAL_DO_INTERSECT_FUNCTION(Circle_2, Triangle_2, 2)

} // namespace CGAL
#endif // CGAL_INTERSECTIONS_2_CIRCLE_2_TRIANGLE_2_H
