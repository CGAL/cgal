// Copyright (c) 2005
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
// Author(s)     : Nico Kruithof

#ifndef CGAL_INTERNAL_INTERSECTIONS_TETRAHEDRON_3_BOUNDED_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_TETRAHEDRON_3_BOUNDED_3_DO_INTERSECT_H

#include <CGAL/Intersections_3/internal/Iso_cuboid_3_Triangle_3_do_intersect.h>
#include <CGAL/Intersections_3/internal/Line_3_Triangle_3_do_intersect.h>
#include <CGAL/Intersections_3/internal/Ray_3_Triangle_3_do_intersect.h>
#include <CGAL/Intersections_3/internal/Segment_3_Triangle_3_do_intersect.h>
#include <CGAL/Intersections_3/internal/Sphere_3_Triangle_3_do_intersect.h>
#include <CGAL/Intersections_3/internal/Tetrahedron_3_Triangle_3_do_intersect.h>
#include <CGAL/Intersections_3/internal/Triangle_3_Triangle_3_do_intersect.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Tetrahedron_3& tet,
             const typename K::Triangle_3& tr,
             const K& k);

// This code is not optimized:
template <class K, class Bounded>
typename K::Boolean
do_intersect_tetrahedron_bounded(const Bounded& tr,
                                 const typename K::Tetrahedron_3& tet,
                                 const typename K::Point_3& p,
                                 const K& k)
{
  typedef typename K::Boolean Boolean;

  CGAL_kernel_precondition(!k.is_degenerate_3_object()(tr));
  CGAL_kernel_precondition(!k.is_degenerate_3_object()(tet));

  Boolean result = false;
  for (int i = 0; i < 4; ++i)
  {
    const Boolean b = do_intersect(tr,
                                   k.construct_triangle_3_object()(tet[i],
                                                                   tet[(i+1)%4],
                                                                   tet[(i+2)%4]),
                                   k);
    if(certainly(b))
      return b;

    if(is_indeterminate(b))
      result = b;
  }

  const Boolean b = k.has_on_bounded_side_3_object()(tet, p);
  if(certainly(b))
    return b;

  if(is_indeterminate(b))
    result = b;

  return result;
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_TETRAHEDRON_3_BOUNDED_3_DO_INTERSECT_H
