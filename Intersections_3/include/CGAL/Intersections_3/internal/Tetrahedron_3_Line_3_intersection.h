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

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_TETRAHEDRON_3_LINE_3_INTERSECTION_H
#define CGAL_INTERNAL_INTERSECTIONS_3_TETRAHEDRON_3_LINE_3_INTERSECTION_H

#include <CGAL/kernel_basic.h>
#include <CGAL/intersections.h>
#include <CGAL/Intersections_3/internal/tetrahedron_lines_intersections_3.h>

namespace CGAL {

namespace Intersections {

namespace internal {

template<class K>
struct Tetrahedron_Line_intersection_3
    :public Tetrahedron_lines_intersection_3_base<K, typename K::Line_3, Tetrahedron_Line_intersection_3<K> >
{
  typedef Tetrahedron_lines_intersection_3_base<K, typename K::Line_3,
  Tetrahedron_Line_intersection_3<K> > Base;
  typedef typename Base::Result_type Result_type;
  Tetrahedron_Line_intersection_3(
      const typename K::Tetrahedron_3& tet,
      const typename K::Line_3& o):Base(tet,o)
  {}

  bool all_inside_test()
  {
    return false;
  }

  bool are_extremities_inside_test()
  {
    return false;
  }

};

//Tetrahedron_3 Line_3
template <class K>
typename Intersection_traits<K, typename K::Tetrahedron_3, typename K::Line_3>::result_type
intersection(
    const typename K::Tetrahedron_3 &tet,
    const typename K::Line_3 &line,
    const K&)
{
  Tetrahedron_Line_intersection_3<K> solver(tet, line);
  solver.do_procede();
  return solver.output;
}

template <class K>
typename Intersection_traits<K, typename K::Tetrahedron_3, typename K::Line_3>::result_type
intersection(
    const typename K::Line_3 &line,
    const typename K::Tetrahedron_3 &tet,
    const K& k)
{
  return intersection(tet, line, k);
}

}}}
#endif // CGAL_INTERNAL_INTERSECTIONS_3_TETRAHEDRON_3_LINE_3_INTERSECTION_H
