// Copyright (c) 2010 GeometryFactory (France).
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

#ifndef CGAL_INTERSECTIONS_3_LINE_3_TETRAHEDRON_3_H
#define CGAL_INTERSECTIONS_3_LINE_3_TETRAHEDRON_3_H

#include <CGAL/Line_3.h>
#include <CGAL/Tetrahedron_3.h>

#include <CGAL/Intersections_3/internal/Tetrahedron_3_Line_3_intersection.h>

namespace CGAL {
  CGAL_DO_INTERSECT_FUNCTION(Line_3, Tetrahedron_3, 3)

  template<typename K>
  typename Intersection_traits<K, typename K::Line_3, typename K::Tetrahedron_3>::result_type
  intersection(const Line_3<K>& line,
               const Tetrahedron_3<K>& tet) {
    return K().intersect_3_object()(line, tet);
  }

  template<typename K>
  typename Intersection_traits<K, typename K::Line_3, typename K::Tetrahedron_3>::result_type
  intersection(const Tetrahedron_3<K>& tet,
      const Line_3<K>& line) {
    return K().intersect_3_object()(tet, line);
  }
}

#endif // CGAL_INTERSECTIONS_3_LINE_3_TETRAHEDRON_3_H
