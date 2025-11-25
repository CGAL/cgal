// Copyright (c) 2021 GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Antonio Gomes, Dmitry Anisimov
//

#ifndef CGAL_BARYCENTRIC_ENUM_3_H
#define CGAL_BARYCENTRIC_ENUM_3_H

#include <CGAL/license/Barycentric_coordinates_3.h>

namespace CGAL {

/*!
  \ingroup PkgBarycentricCoordinates3Ref

  The namespace `Barycentric_coordinates` contains implementations of all
  generalized barycentric coordinates: 2D, 3D, related enumerations, etc.
*/
namespace Barycentric_coordinates {

/// \name Computation Policies
/// @{

/*!
  `Computation_policy_3` provides a way to choose an asymptotic time complexity
  of the algorithm and its precision for computing 3D barycentric weights and coordinates.
*/
enum class Computation_policy_3 {

  /*!
    Computation has a linear time complexity with respect to the number of
    face vertices, but may suffer imprecisions near the face boundary.
    No extra checks are carried out.
  */
  FAST = 0,

    /*!
    Computation has a linear time complexity with respect to the number of
    face vertices, but may suffer imprecisions near the face boundary. In
    addition, we check a position of the query point with respect to the face
    and use different computation strategies for different positions.
  */
  FAST_WITH_EDGE_CASES = 1
};

/// @}

} // namespace Barycentric_coordinates
} // namespace CGAL

#endif // CGAL_BARYCENTRIC_ENUM_3_H
