// Copyright (c) 2018 GeometryFactory (France).
// All rights reserved.
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
// Author(s)     : Andreas Fabri


#ifndef CGAL_INTERNAL_INTERSECTIONS_2_BBOX_2_CIRCLE_2_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_2_BBOX_2_CIRCLE_2_DO_INTERSECT_H

#include <CGAL/Circle_2.h>
#include <CGAL/Bbox_2.h>

#include <CGAL/number_utils.h>


namespace CGAL {
namespace Intersections {
namespace internal {

// Circle_2 is not a disk, thus if the box is contained within the circle, there is no intersection.
template <class K, class Box2>
bool do_intersect_circle_box_2(const typename K::Circle_2& circle,
                               const Box2& bbox,
                               const K&)
{
  typedef typename K::FT                                          FT;
  typedef typename K::Point_2                                     Point;

  Point center = circle.center();

  // Check that the minimum distance to the box is smaller than the radius, otherwise there is
  // no intersection.
  FT distance = FT(0);
  if(center.x() < FT(bbox.xmin()))
  {
    FT d = FT(bbox.xmin()) - center.x();
    distance += d * d;
  }
  else if(center.x() > FT(bbox.xmax()))
  {
    FT d = center.x() - FT(bbox.xmax());
    distance += d * d;
  }

  if(center.y() < FT(bbox.ymin()))
  {
    FT d = FT(bbox.ymin()) - center.y();
    distance += d * d;
  }
  else if(center.y() > FT(bbox.ymax()))
  {
    FT d = center.y() - FT(bbox.ymax());
    distance += d * d;
  }

  // Note that with the way the distance above is computed, the distance is '0' if the box strictly
  // contains the circle. But since we use '>', we don't exit
  if(distance > circle.squared_radius())
    return false;

  // Check that the maximum distance between the center of the circle and the box is not (strictly)
  // smaller than the radius of the center, otherwise the box is entirely contained.
  distance = FT(0);
  if(center.x() <= FT(bbox.xmin() + bbox.xmax()) / FT(2))
  {
    FT d = FT(bbox.xmax()) - center.x();
    distance += d * d;
  }
  else
  {
    FT d = center.x() - FT(bbox.xmin());
    distance += d * d;
  }

  if(center.y() < FT(bbox.ymin() + bbox.ymax()) / FT(2))
  {
    FT d = FT(bbox.ymax()) - center.y();
    distance += d * d;
  }
  else
  {
    FT d = center.y() - FT(bbox.ymin());
    distance += d * d;
  }

  return (distance >= circle.squared_radius());
}

template <class K>
bool do_intersect(const CGAL::Bbox_2& bbox,
                  const typename K::Circle_2& circle,
                  const K&)
{
  return do_intersect_circle_box_2(circle, bbox, K());
}


template <class K>
bool do_intersect(const typename K::Circle_2& circle,
                  const CGAL::Bbox_2& bbox,
                  const K&)
{
  return do_intersect_circle_box_2(circle, bbox, K());
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif  // CGAL_INTERNAL_INTERSECTIONS_2_BBOX_2_CIRCLE_2_DO_INTERSECT_H
