// Copyright (c) 2019 GeometryFactory(France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Maxime Gimeno
//

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_TETRAHEDRON_3_RAY_3_INTERSECTION_H
#define CGAL_INTERNAL_INTERSECTIONS_3_TETRAHEDRON_3_RAY_3_INTERSECTION_H

#include <CGAL/kernel_basic.h>
#include <CGAL/intersections.h>
#include <CGAL/Intersections_3/internal/tetrahedron_lines_intersections_3.h>

namespace CGAL {

namespace Intersections {

namespace internal {

template<class K>
struct Tetrahedron_ray_intersection_3
    :public Tetrahedron_lines_intersection_3_base<K,
      typename K::Ray_3, Tetrahedron_ray_intersection_3<K> >
{
  typedef typename K::Ray_3 O;
  typedef Tetrahedron_lines_intersection_3_base<K, typename K::Ray_3,
    Tetrahedron_ray_intersection_3<K> > Base;
  typedef typename Base::Result_type Result_type;

  Tetrahedron_ray_intersection_3(const typename K::Tetrahedron_3& tet,
                                 const O& o):Base(tet,o) {}

  bool are_extremities_inside_test()
  {
    //If one extremity is inside tet : return a segment
    if(this->tet.has_on_bounded_side(this->o.source()))
    {
      typename K::Segment_3 result(this->o.source(), this->res_points.front());
      this->output = Result_type(std::forward<typename K::Segment_3>(result));
      return true;
    }
    return false;
  }
};


//Tetrahedron_3 Ray_3
template <class K>
typename Intersection_traits<K, typename K::Tetrahedron_3, typename K::Ray_3>::result_type
intersection(
    const typename K::Tetrahedron_3 &tet,
    const typename K::Ray_3 &ray,
    const K&)
{
  Tetrahedron_ray_intersection_3<K> solver(tet, ray);
  solver.do_procede();
  return solver.output;
}

template <class K>
typename Intersection_traits<K, typename K::Tetrahedron_3, typename K::Ray_3>::result_type
intersection(
    const typename K::Ray_3 &seg,
    const typename K::Tetrahedron_3 &tet,
    const K& k)
{
  return intersection(tet, seg, k);
}

}}}
#endif // CGAL_INTERNAL_INTERSECTIONS_3_TETRAHEDRON_3_RAY_3_INTERSECTION_H
