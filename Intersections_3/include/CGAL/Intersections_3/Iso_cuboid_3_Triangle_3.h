// Copyright (c) 1997-2021
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).
// GeometryFactory (France)
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

#ifndef CGAL_INTERSECTIONS_3_ISO_CUBOID_3_TRIANGLE_3_H
#define CGAL_INTERSECTIONS_3_ISO_CUBOID_3_TRIANGLE_3_H

#include <CGAL/Intersection_traits.h>
#include <CGAL/Intersections_3/internal/Iso_cuboid_3_Triangle_3_do_intersect.h>
#include <CGAL/Intersections_3/internal/Iso_cuboid_3_Triangle_3_intersection.h>

#include <CGAL/Iso_cuboid_3.h>
#include <CGAL/Triangle_3.h>

namespace CGAL {

CGAL_DO_INTERSECT_FUNCTION(Iso_cuboid_3, Triangle_3, 3)
CGAL_INTERSECTION_FUNCTION(Iso_cuboid_3, Triangle_3, 3)

} // namespace CGAL

#endif // CGAL_INTERSECTIONS_3_BBOX_3_TRIANGLE_3_H
