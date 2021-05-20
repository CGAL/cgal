// Copyright (c) 2019 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Dmitry Anisimov, David Bommes, Kai Hormann, Pierre Alliez
//

#ifndef CGAL_BARYCENTRIC_ENUM_2_H
#define CGAL_BARYCENTRIC_ENUM_2_H

#include <CGAL/license/Barycentric_coordinates_2.h>

namespace CGAL {

/*!
  \ingroup PkgBarycentricCoordinates2Ref

  The namespace `Barycentric_coordinates` contains implementations of all
  generalized barycentric coordinates: 2D, 3D, related enumerations, etc.
*/
namespace Barycentric_coordinates {

/// \name Computation Policies
/// @{

/*!
  `Computation_policy_2` provides a way to choose an asymptotic time complexity
  of the algorithm and its precision for computing 2D barycentric weights and coordinates.
*/
enum class Computation_policy_2 {

  /*!
    Computation is very precise but has a quadratic time complexity with respect
    to the number of polygon vertices. In addition, we check a position of
    the query point with respect to the polygon and use different computation
    strategies for different positions. This is the default policy.
  */
  PRECISE_WITH_EDGE_CASES = 0,

  /*!
    Computation is very precise but has a quadratic time complexity with respect
    to the number of polygon vertices. No extra checks are carried out.
  */
  PRECISE = 1,

  /*!
    Computation has a linear time complexity with respect to the number of
    polygon vertices, but may suffer imprecisions near the polygon boundary. In
    addition, we check a position of the query point with respect to the polygon
    and use different computation strategies for different positions.
  */
  FAST_WITH_EDGE_CASES = 2,

  /*!
    Computation has a linear time complexity with respect to the number of
    polygon vertices, but may suffer imprecisions near the polygon boundary.
    No extra checks are carried out.
  */
  FAST = 3
};

/// @}

} // namespace Barycentric_coordinates
} // namespace CGAL

#endif // CGAL_BARYCENTRIC_ENUM_2_H
