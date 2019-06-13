// Copyright (c) 2019 GeometryFactory(France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
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



//Tetrahedron_3 Line_3
template <class K>
typename Intersection_traits<K, typename K::Tetrahedron_3, typename K::Line_3>::result_type
intersection(
    const typename K::Tetrahedron_3 &tet,
    const typename K::Line_3 &line,
    const K&)
{
  Tetrahedron_lines_intersection_3_base<K, typename K::Line_3> solver(tet, line);
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
