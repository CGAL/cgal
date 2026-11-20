// Copyright (c) 2026 GeometryFactory Sarl
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

#ifndef CGAL_INTERSECTIONS_3_BALL_3_ISO_CUBOID_3_H
#define CGAL_INTERSECTIONS_3_BALL_3_ISO_CUBOID_3_H

#include <CGAL/Intersection_traits_3.h>
#include <CGAL/Intersections_3/internal/Ball_3_Iso_cuboid_3_do_intersect.h>

#include <CGAL/Ball_3.h>
#include <CGAL/Iso_cuboid_3.h>

namespace CGAL {

CGAL_DO_INTERSECT_FUNCTION(Ball_3, Iso_cuboid_3, 3)

} // namespace CGAL

#endif // CGAL_INTERSECTIONS_3_BALL_3_ISO_CUBOID_3_H
