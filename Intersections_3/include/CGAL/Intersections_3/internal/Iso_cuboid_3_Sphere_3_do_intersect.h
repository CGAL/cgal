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

#include <CGAL/Coercion_traits.h>
#include <CGAL/number_utils.h>
#include <CGAL/Uncertain.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K, class BFT> // Iso_cuboid_3 or Bbox_3
bool do_intersect_sphere_box_3(const typename K::Sphere_3& sphere,
                               const BFT bxmin, const BFT bymin, const BFT bzmin,
                               const BFT bxmax, const BFT bymax, const BFT bzmax,
                               const K&)
{
  typedef typename K::FT SFT;
  typedef typename Coercion_traits<SFT, BFT>::Type FT;
  typedef typename K::Point_3 Point;

  typename Coercion_traits<SFT, BFT>::Cast to_FT;

  FT d = FT(0);
  FT distance = FT(0);
  FT sr = sphere.squared_radius();

  const Point& center = sphere.center();

  if(compare(center.x(), bxmin) == SMALLER)
  {
    d = to_FT(bxmin) - to_FT(center.x());
    d = square(d);
    if(certainly(d > sr))
      return false;

    distance = d;
  }
  else if(compare(center.x(), bxmax) == LARGER)
  {
    d = to_FT(center.x()) - to_FT(bxmax);
    d = square(d);
    if(certainly(d > sr))
      return false;

    distance = d;
  }

  if(compare(center.y(), bymin) == SMALLER)
  {
    d = to_FT(bymin) - to_FT(center.y());
    d = square(d);
    if(certainly(d > sr))
      return false;

    distance += d;
  }
  else if(compare(center.y(), bymax) == LARGER)
  {
    d = to_FT(center.y()) - to_FT(bymax);
    d = square(d);
    if(certainly(d > sr))
      return false;

    distance += d;
  }

  if(compare(center.z(), bzmin) == SMALLER)
  {
    d = to_FT(bzmin) - to_FT(center.z());
    d = square(d);
    distance += d;
  }
  else if(compare(center.z(), bzmax) == LARGER)
  {
    d = to_FT(center.z()) - to_FT(bzmax);
    d = square(d);
    distance += d;
  }

  return (distance <= sr);
}

template <class K>
bool do_intersect(const typename K::Sphere_3& sphere,
                  const typename K::Iso_cuboid_3& ic,
                  const K& k)
{
  return do_intersect_sphere_box_3(sphere,
                                   (ic.min)().x(), (ic.min)().y(), (ic.min)().z(),
                                   (ic.max)().x(), (ic.max)().y(), (ic.max)().z(),
                                   k);
}

template <class K>
bool do_intersect(const typename K::Iso_cuboid_3& ic,
                  const typename K::Sphere_3& sphere,
                  const K& k)
{
  return do_intersect(sphere, ic, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_ISO_CUBOID_3_SPHERE_3_DO_INTERSECT_H
