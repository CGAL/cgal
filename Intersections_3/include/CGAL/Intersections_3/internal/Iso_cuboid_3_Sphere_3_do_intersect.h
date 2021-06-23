// Copyright (c) 2008  INRIA Sophia-Antipolis (France), ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Camille Wormser, Jane Tournois, Pierre Alliez

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_ISO_CUBOID_3_SPHERE_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_3_ISO_CUBOID_3_SPHERE_3_DO_INTERSECT_H

#include <CGAL/number_utils.h>
#include <CGAL/Uncertain.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K, class Box3> // Iso_cuboid_3 or Bbox_3
bool do_intersect_sphere_box_3(const typename K::Sphere_3& sphere,
                               const Box3& bbox,
                               const K&)
{
  typedef typename K::FT FT;
  typedef typename K::Point_3 Point;

  FT d = FT(0);
  FT distance = FT(0);
  FT sr = sphere.squared_radius();

  const Point& center = sphere.center();

  if(center.x() < FT{bbox.xmin()})
  {
    d = FT{bbox.xmin()} - center.x();
    d = square(d);
    if(certainly(d > sr))
      return false;

    distance = d;
  }
  else if(center.x() > FT{bbox.xmax()})
  {
    d = center.x() - FT{bbox.xmax()};
    d = square(d);
    if(certainly(d > sr))
      return false;

    distance = d;
  }

  if(center.y() < FT{bbox.ymin()})
  {
    d = FT{bbox.ymin()} - center.y();
    d = square(d);
    if(certainly(d > sr))
      return false;

    distance += d;
  }
  else if(center.y() > FT{bbox.ymax()})
  {
    d = center.y() - FT{bbox.ymax()};
    d = square(d);
    if(certainly(d > sr))
      return false;

    distance += d;
  }

  if(center.z() < FT{bbox.zmin()})
  {
    d = FT{bbox.zmin()} - center.z();
    d = square(d);
    distance += d;
  }
  else if(center.z() > FT{bbox.zmax()})
  {
    d = center.z() - FT{bbox.zmax()};
    d = square(d);
    distance += d;
  }

//  For unknown reasons this causes a syntax error on VC2005
//  but compiles fine on Linux and MAC
//
//  int i;
//  for(i = 0; i < 3; ++i)
//  {
//    if(center[i] < (FT)bbox.min(i))
//    {
//      d = (FT)bbox.min(i) - center[i];
//      distance += d * d;
//    }
//    else if(center[i] > (FT)bbox.max(i))
//    {
//      d = center[i] - (FT)bbox.max(i);
//      distance += d * d;
//    }
//  }

  return distance <= sr;
}

template <class K>
bool do_intersect(const typename K::Iso_cuboid_3& ic,
                  const typename K::Sphere_3& sphere,
                  const K& k)
{
  return do_intersect_sphere_box_3(sphere, ic, k);
}

template <class K>
bool do_intersect(const typename K::Sphere_3& sphere,
                  const typename K::Iso_cuboid_3& ic,
                  const K& k)
{
  return do_intersect_sphere_box_3(sphere, ic, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_ISO_CUBOID_3_SPHERE_3_DO_INTERSECT_H
