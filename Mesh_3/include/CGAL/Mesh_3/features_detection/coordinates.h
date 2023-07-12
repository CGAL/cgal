// Copyright (c) 2022 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sebastien Loriot
//
//******************************************************************************
//
//******************************************************************************


#ifndef CGAL_MESH_3_FEATURES_DETECTION_COORDINATES_H
#define CGAL_MESH_3_FEATURES_DETECTION_COORDINATES_H

#include <CGAL/license/Mesh_3.h>

#include <array>

namespace CGAL
{
namespace Mesh_3
{
namespace internal
{
  using Coordinates = std::array<int, 3>;
  constexpr Coordinates coordinates[8] = { {0, 0, 0},
                                           {1, 0, 0},
                                           {0, 1, 0},
                                           {1, 1, 0},
                                           {0, 0, 1},
                                           {1, 0, 1},
                                           {0, 1, 1},
                                           {1, 1, 1} };

  inline Coordinates minus(const Coordinates& b, const Coordinates& a) {
    return { b[0] - a[0], b[1] - a[1], b[2] - a[2] };
  }

  inline Coordinates cross(Coordinates a, Coordinates b) {
    return { a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0] };
  }

  inline Coordinates square(Coordinates c) {
    return { c[0] * c[0], c[1] * c[1], c[2] * c[2] };
  }

  inline int dist(Coordinates a, Coordinates b) {
    auto s = square(minus(b, a));
    return s[0] + s[1] + s[2];
  }

}//end namespace internal
}//end namespace Mesh_3
}//end namespace CGAL

#endif // CGAL_MESH_3_FEATURES_DETECTION_COORDINATES_H
