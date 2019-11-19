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

#ifndef CGAL_INTERSECTIONS_3_ISO_CUBOID_3_RAY_3_H
#define CGAL_INTERSECTIONS_3_ISO_CUBOID_3_RAY_3_H

#include <CGAL/Iso_cuboid_3.h>
#include <CGAL/Ray_3.h>

#include <CGAL/Intersections_3/internal/intersection_3_1_impl.h>
#include <CGAL/Intersections_3/internal/Iso_cuboid_3_Ray_3_do_intersect.h>

namespace CGAL {
  CGAL_DO_INTERSECT_FUNCTION(Iso_cuboid_3,Ray_3, 3)
  CGAL_INTERSECTION_FUNCTION(Ray_3, Iso_cuboid_3, 3)
}

#endif // CGAL_INTERSECTIONS_3_ISSO_CUBOID_3_RAY_3_H
