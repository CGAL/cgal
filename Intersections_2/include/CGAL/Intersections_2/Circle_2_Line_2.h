// Copyright (c) 2000
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Geert-Jan Giezeman


#ifndef CGAL_INTERSECTIONS_2_CIRCLE_2_LINE_2_H
#define CGAL_INTERSECTIONS_2_CIRCLE_2_LINE_2_H

#include <CGAL/Circle_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/squared_distance_2_1.h>
#include <CGAL/Intersection_traits_2.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
bool
do_intersect(const typename K::Circle_2 & c,
             const typename K::Line_2& l,
             const K&)
{
  return squared_distance(c.center(), l) <= c.squared_radius();
}

template <class K>
bool
do_intersect(const typename K::Line_2& l,
             const typename K::Circle_2 & c,
             const K&)
{
  return squared_distance(c.center(), l) <= c.squared_radius();
}

} // namespace internal
} // namespace Intersections

CGAL_DO_INTERSECT_FUNCTION(Circle_2, Line_2, 2)

} // namespace CGAL

#endif
