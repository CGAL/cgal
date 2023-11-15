// Copyright (c) 2019 GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Simon Giraudot, Dmitry Anisimov

#ifndef CGAL_KSR_ENUM_H
#define CGAL_KSR_ENUM_H

#include <CGAL/license/Kinetic_shape_partition.h>

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

  enum class Visibility_label {

    /// Outside the object.
    OUTSIDE = 0,

    /// Inside the object.
    INSIDE = 1
  };

} // namespace KSR
} // namespace CGAL

#endif // CGAL_KSR_ENUM_H
