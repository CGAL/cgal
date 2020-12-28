// Copyright (c) 2019 GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_KSR_ENUM_H
#define CGAL_KSR_ENUM_H

// #include <CGAL/license/Kinetic_shape_reconstruction.h>

namespace CGAL {
namespace KSR {

  enum class Semantic_label {

    // Ground points.
		GROUND = 0,

    // Items treated as vegetation.
    VEGETATION = 1,

    // Items treated as building boundary, e.g. walls.
		BUILDING_BOUNDARY = 2,

    // Items treated as building interior, e.g. roofs.
    BUILDING_INTERIOR = 3,

    // Any item that is not handled by the algorithm.
    UNCLASSIFIED = 4

	};

  enum class Planar_shape_type {

    // Convex hull shape type.
    CONVEX_HULL = 0,

    // Rectangle shape type.
    RECTANGLE = 1
  };

} // namespace KSR
} // namespace CGAL

#endif // CGAL_KSR_ENUM_H
