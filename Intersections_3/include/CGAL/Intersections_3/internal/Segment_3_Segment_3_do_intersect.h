// Copyright (c) 2003  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Philippe Guigue

#ifndef CGAL_INTERNAL_INTERSECTIONS_SEGMENT_3_SEGMENT_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_SEGMENT_3_SEGMENT_3_DO_INTERSECT_H

#include <CGAL/Intersections_3/internal/Line_3_Line_3_do_intersect.h>

#include <CGAL/enum.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
inline
bool
do_intersect(const typename K::Segment_3& s1,
             const typename K::Segment_3& s2,
             const K& k)
{
  CGAL_precondition(!s1.is_degenerate() && !s2.is_degenerate());

  bool b = do_intersect(s1.supporting_line(), s2.supporting_line(), k);
  if(b)
  {
    // supporting_line intersects: points are coplanar
    typename K::Coplanar_orientation_3 cpl_orient=k.coplanar_orientation_3_object();
    ::CGAL::Orientation or1 = cpl_orient(s1[0], s1[1], s2[0]);
    ::CGAL::Orientation or2 = cpl_orient(s1[0], s1[1], s2[1]);

    if(or1 == COLLINEAR && or2 == COLLINEAR)
    {
      // segments are collinear
      typename K::Collinear_are_ordered_along_line_3 cln_order = k.collinear_are_ordered_along_line_3_object();
      return (cln_order(s1[0], s2[0], s1[1]) ||
              cln_order(s1[0], s2[1], s1[1]) ||
              cln_order(s2[0], s1[0], s2[1]));
    }

    if(or1 != or2)
    {
      or1 = cpl_orient(s2[0], s2[1], s1[0]);
      return (or1 == COLLINEAR || or1 != cpl_orient(s2[0], s2[1], s1[1]));
    }
  }

  return false;
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif //CGAL_INTERNAL_INTERSECTIONS_SEGMENT_3_SEGMENT_3_DO_INTERSECT_H
