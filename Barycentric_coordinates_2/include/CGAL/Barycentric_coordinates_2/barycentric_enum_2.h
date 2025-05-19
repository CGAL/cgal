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

namespace CGAL {
namespace Barycentric_coordinates {

#if !defined(CGAL_NO_DEPRECATED_CODE) || defined(DOXYGEN_RUNNING)

/// \name Locations of a Query Point
/// @{

/// Query_point_location is enumeration with possible locations of a query point provided by the user.
/// \deprecated This part of the package is deprecated since the version 5.4 of \cgal.
enum
#ifndef DOXYGEN_RUNNING
CGAL_DEPRECATED_MSG("This part of the package is deprecated since the version 5.4 of CGAL!")
#endif
Query_point_location
{
  /// Location is not known apriori and is defined automatically by the algorithm.
  UNSPECIFIED_LOCATION,

  /// Query point is located at one of the polygon's vertices.
  ON_VERTEX,

  /// Query point is located on the boundary of the polygon.
  ON_BOUNDARY,

  /// Query point is located inside the polygon, excluding the boundary.
  ON_BOUNDED_SIDE,

  /// Query point is located outside the polygon, excluding the boundary.
  ON_UNBOUNDED_SIDE
};

/// @}

/// \name Types of an Algorithm
/// @{

/// Type_of_algorithm is enumeration with possible algorithms to compute coordinates.
/// \deprecated This part of the package is deprecated since the version 5.4 of \cgal.
enum
#ifndef DOXYGEN_RUNNING
CGAL_DEPRECATED_MSG("This part of the package is deprecated since the version 5.4 of CGAL!")
#endif
Type_of_algorithm
{
  /// A default slow algorithm, which is as precise as possible.
  PRECISE,

  /// A fast algorithm, which is less precise but much faster.
  FAST
};

/// @}

/// \cond SKIP_IN_MANUAL

// Types of a Polygon

// Type_of_polygon is enumeration with possible types of the input polygon.
// It is used internally to precondition coordinates.
enum Type_of_polygon {

  // Concave polygon = non-convex polygon.
  CONCAVE,

  // This is a convex polygon with collinear vertices.
  WEAKLY_CONVEX,

  // This is a convex polygon without collinear vertices.
  STRICTLY_CONVEX
};

/// \endcond

#endif // CGAL_NO_DEPRECATED_CODE

} // namespace Barycentric_coordinates
} // namespace CGAL

#endif // CGAL_BARYCENTRIC_ENUM_2_H
