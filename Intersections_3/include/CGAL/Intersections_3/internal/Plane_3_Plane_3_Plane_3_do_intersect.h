// Copyright (c) 2008 INRIA(France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Geert-Jan Giezeman,
//                 Andreas Fabri

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_PLANE_3_PLANE_3_PLANE_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_3_PLANE_3_PLANE_3_PLANE_3_DO_INTERSECT_H

#include <CGAL/determinant.h>
#include <CGAL/number_utils.h>
#include <CGAL/rank.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
inline bool
do_intersect(const typename K::Plane_3& plane1,
             const typename K::Plane_3& plane2,
             const typename K::Plane_3& plane3,
             const K&)
{
  typedef typename K::RT RT;

  if(!is_zero(determinant(plane1.a(), plane1.b(), plane1.c(),
                          plane2.a(), plane2.b(), plane2.c(),
                          plane3.a(), plane3.b(), plane3.c())))
    return true;

  int pcount = 0;
  bool b12, b13,b23;
  if((b12 = parallel(plane1, plane2))) pcount++;
  if((b13 = parallel(plane1, plane3))) pcount++;
  if((b23 = parallel(plane2, plane3))) pcount++;

  if(pcount == 3)
  {
    return (((plane1 == plane2) || (plane1 == plane2.opposite())) &&
            ((plane1 == plane3) || (plane1 == plane3.opposite())));
  }

  if(pcount == 1)
  {
    if(b12 && ((plane1 == plane2)||(plane1 == plane2.opposite()))) return true;
    if(b13 && ((plane1 == plane3)||(plane1 == plane3.opposite()))) return true;
    if(b23 && ((plane2 == plane3)||(plane2 == plane3.opposite()))) return true;
  }

  int rd = rank_34<RT>(plane1.a(), plane1.b(), plane1.c(), plane1.d(),
                       plane2.a(), plane2.b(), plane2.c(), plane2.d(),
                       plane3.a(), plane3.b(), plane3.c(), plane3.d());

  return (rd == 2);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_PLANE_3_PLANE_3_PLANE_3_DO_INTERSECT_H
