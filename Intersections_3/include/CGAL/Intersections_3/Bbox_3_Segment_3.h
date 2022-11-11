// Copyright (c) 1997-2021
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).
// GeometryFactory (France)
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

#ifndef CGAL_INTERSECTIONS_3_BBOX_3_SEGMENT_3_H
#define CGAL_INTERSECTIONS_3_BBOX_3_SEGMENT_3_H

#include <CGAL/Intersection_traits_3.h>
#include <CGAL/Intersections_3/internal/Bbox_3_Segment_3_do_intersect.h>
#include <CGAL/Intersections_3/internal/Bbox_3_Segment_3_intersection.h>

#include <CGAL/Bbox_3.h>
#include <CGAL/Segment_3.h>

namespace CGAL {

template<typename K>
bool do_intersect(const CGAL::Bbox_3& box,
                  const Segment_3<K>& s)
{
  return K().do_intersect_3_object()(box, s);
}

template<typename K>
bool do_intersect(const Segment_3<K>& s,
                  const CGAL::Bbox_3& box)
{
  return K().do_intersect_3_object()(s, box);
}

template<typename K>
typename Intersection_traits<K, typename K::Segment_3, Bbox_3>::result_type
intersection(const CGAL::Bbox_3& box,
             const Segment_3<K>& s)
{
  return K().intersect_3_object()(box, s);
}

template<typename K>
typename Intersection_traits<K, typename K::Segment_3, Bbox_3>::result_type
intersection(const Segment_3<K>& s,
             const CGAL::Bbox_3& box)
{
  return K().intersect_3_object()(s, box);
}

} // namespace CGAL

#endif // CGAL_INTERSECTIONS_3_BBOX_3_SEGMENT_3_H
