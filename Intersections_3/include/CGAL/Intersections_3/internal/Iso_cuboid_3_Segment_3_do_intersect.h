// Copyright (c) 2008 ETH Zurich (Switzerland)
// Copyright (c) 2008-2009 INRIA Sophia-Antipolis (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Camille Wormser, Jane Tournois, Pierre Alliez, Stephane Tayeb

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_ISO_CUBOID_3_SEGMENT_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_3_ISO_CUBOID_3_SEGMENT_3_DO_INTERSECT_H

// inspired from http://cag.csail.mit.edu/~amy/papers/box-jgt.pdf

#include <CGAL/Intersections_3/internal/Bbox_3_Segment_3_do_intersect.h>
// for CGAL::internal::do_intersect_bbox_segment_aux

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
bool do_intersect(const typename K::Segment_3& seg,
                  const typename K::Iso_cuboid_3& ic,
                  const K&)
{
  typedef typename K::FT FT;
  typedef typename K::Point_3 Point_3;

  const Point_3& source = seg.source();
  const Point_3& target = seg.target();

  return do_intersect_bbox_segment_aux
    <FT,FT,
     true,  // bounded at t=0
     true,  // bounded at t=1
     false> // do not use static filters
    (
     source.x(), source.y(), source.z(),
     target.x(), target.y(), target.z(),
     (ic.min)().x(), (ic.min)().y(), (ic.min)().z(),
     (ic.max)().x(), (ic.max)().y(), (ic.max)().z()
     );
}

template <class K>
bool do_intersect(const typename K::Iso_cuboid_3& ic,
                  const typename K::Segment_3& seg,
                  const K& k)
{
  return do_intersect(seg, ic, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_ISO_CUBOID_3_SEGMENT_3_DO_INTERSECT_H
