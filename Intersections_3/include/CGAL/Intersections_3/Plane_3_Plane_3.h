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

#ifndef CGAL_INTERSECTIONS_3_PLANE_3_PLANE_3_H
#define CGAL_INTERSECTIONS_3_PLANE_3_PLANE_3_H

#include <CGAL/Plane_3.h>

#include <CGAL/Intersections_3/internal/intersection_3_1_impl.h>

namespace CGAL {
CGAL_INTERSECTION_FUNCTION_SELF(Plane_3, 3)
CGAL_DO_INTERSECT_FUNCTION_SELF(Plane_3, 3)


template < class K >
inline
boost::optional<typename K::Point_3>
intersection_point_for_polyhedral_envelope(const Plane_3<K>& p0, const Plane_3<K>& p1, const Plane_3<K>& p2)
{
  return K().intersect_point_3_for_polyhedral_envelope_object()(p0, p1, p2);
}

}

#endif // CGAL_INTERSECTIONS_3_PLANE_3_PLANE_3_H
