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

#ifndef CGAL_INTERSECTIONS_3_LINE_PLANE_3_H
#define CGAL_INTERSECTIONS_3_LINE_PLANE_3_H

#include <CGAL/Line_3.h>
#include <CGAL/Plane_3.h>

#include <CGAL/Intersections_3/internal/intersection_3_1_impl.h>

namespace CGAL {
CGAL_INTERSECTION_FUNCTION(Plane_3, Line_3, 3)
CGAL_DO_INTERSECT_FUNCTION(Plane_3, Line_3, 3)
}

#endif // CGAL_INTERSECTIONS_3_LINE_PLANE_3_H
