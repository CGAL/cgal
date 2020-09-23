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

#ifndef CGAL_TRIANGLE_3_PLANE_3_DO_INTERSECT_H
#define CGAL_TRIANGLE_3_PLANE_3_DO_INTERSECT_H

#include <CGAL/enum.h>
#include <CGAL/kernel_assertions.h>
#include <CGAL/Intersection_traits_3.h>

namespace CGAL {

  template <class K>
  class Triangle_3;

  template <class K>
  class Plane_3;

namespace Intersections {

namespace internal {

template <class K>
bool do_intersect(const typename K::Triangle_3 &t,
                  const typename K::Plane_3   &h,
                  const K & k)
{

  CGAL_kernel_precondition( ! k.is_degenerate_3_object()(t)) ;
  CGAL_kernel_precondition( ! k.is_degenerate_3_object()(h)) ;


  typename K::Construct_vertex_3 vertex_on =
    k.construct_vertex_3_object();

  typename K::Oriented_side_3 oriented_side =
    k.oriented_side_3_object();



  switch ( oriented_side(h,vertex_on(t,0)) ) {
  case ON_POSITIVE_SIDE:
    return oriented_side(h,vertex_on(t,1)) != ON_POSITIVE_SIDE
      || oriented_side(h,vertex_on(t,2)) != ON_POSITIVE_SIDE;
  case ON_NEGATIVE_SIDE:
    return oriented_side(h,vertex_on(t,1)) != ON_NEGATIVE_SIDE
      || oriented_side(h,vertex_on(t,2)) != ON_NEGATIVE_SIDE ;
  case ON_ORIENTED_BOUNDARY:
    return true;
  default:// should not happen.
    CGAL_kernel_assertion(false);
    return false;
  }
}


template <class K>
inline
bool do_intersect(const typename K::Plane_3   &h,
                  const typename K::Triangle_3 &t,
                  const K & k)
{
  return do_intersect(t, h, k);
}


} // namespace internal
} // namespace Intersections
} //namespace CGAL

#endif //CGAL_TRIANGLE_3_PLANE_3_DO_INTERSECT_H
