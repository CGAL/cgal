// Copyright (c) 2018 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri
//

#ifndef CGAL_INTERSECTIONS_2_ISO_RECTANGLE_2_CIRCLE_2_H
#define CGAL_INTERSECTIONS_2_ISO_RECTANGLE_2_CIRCLE_2_H

#include <CGAL/Circle_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Intersection_traits_2.h>

namespace CGAL {
namespace Intersections {
namespace internal {

// Circle_2 is not a disk, thus if the box is contained within the circle, there is no intersection.
template <class K>
bool do_intersect_circle_iso_rectangle_2(const typename K::Circle_2& circle,
                                         const typename K::Iso_rectangle_2& rec,
                                         const K&)
{
  typedef typename K::FT                                          FT;
  typedef typename K::Point_2                                     Point;

  Point center = circle.center();

  // Check that the minimum distance to the box is smaller than the radius, otherwise there is
  // no intersection.
  FT distance = FT(0);
  if(center.x() < rec.xmin())
  {
    FT d = rec.xmin() - center.x();
    distance += d * d;
  }
  else if(center.x() > rec.xmax())
  {
    FT d = center.x() - rec.xmax();
    distance += d * d;
  }

  if(center.y() < rec.ymin())
  {
    FT d = rec.ymin() - center.y();
    distance += d * d;
  }
  else if(center.y() > rec.ymax())
  {
    FT d = center.y() - rec.ymax();
    distance += d * d;
  }

  // Note that with the way the distance above is computed, the distance is '0' if the box strictly
  // contains the circle. But since we use '>', we don't exit
  if(distance > circle.squared_radius())
    return false;

  // Check that the maximum distance between the center of the circle and the box is not (strictly)
  // smaller than the radius of the center, otherwise the box is entirely contained.
  distance = FT(0);
  if(center.x() <= (rec.xmin() + rec.xmax()) / FT(2))
  {
    FT d = rec.xmax() - center.x();
    distance += d * d;
  }
  else
  {
    FT d = center.x() - rec.xmin();
    distance += d * d;
  }

  if(center.y() < (rec.ymin() + rec.ymax()) / FT(2))
  {
    FT d = rec.ymax() - center.y();
    distance += d * d;
  }
  else
  {
    FT d = center.y() - rec.ymin();
    distance += d * d;
  }

  return (distance >= circle.squared_radius());
}

template <class K>
bool do_intersect(const typename K::Iso_rectangle_2& rec,
                  const typename K::Circle_2& circle,
                  const K&)
{
  return do_intersect_circle_iso_rectangle_2(circle, rec, K());
}


template <class K>
bool do_intersect(const typename K::Circle_2& circle,
                  const typename K::Iso_rectangle_2& rec,
                  const K&)
{
  return do_intersect_circle_iso_rectangle_2(circle, rec, K());
}

} // namespace internal
} // namespace Intersections

CGAL_DO_INTERSECT_FUNCTION(Iso_rectangle_2, Circle_2, 2)

}

#endif // CGAL_INTERSECTIONS_2_ISO_RECTANGLE_2_CIRCLE_2_H
