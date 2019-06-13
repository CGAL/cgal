// Copyright (c) 2018 GeometryFactory (France).
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
// Author(s)     : Andreas Fabri
//

#ifndef CGAL_INTERSECTIONS_2_BBOX_2_CIRCLE_2_H
#define CGAL_INTERSECTIONS_2_BBOX_2_CIRCLE_2_H

#include <CGAL/Bbox_2.h>
#include <CGAL/Circle_2.h>

#include <CGAL/Intersections_2/internal/Bbox_2_Circle_2_do_intersect.h>

namespace CGAL {

template<typename K>
bool do_intersect(const CGAL::Bbox_2& a,
                  const Circle_2<K>& b) {
  return K().do_intersect_2_object()(a, b);
}

template<typename K>
bool do_intersect(const Circle_2<K>& a,
                  const CGAL::Bbox_2& b) {
  return K().do_intersect_2_object()(a, b);
}


}
#endif // CGAL_INTERSECTIONS_2_BBOX_2_CIRCLE_2_H
