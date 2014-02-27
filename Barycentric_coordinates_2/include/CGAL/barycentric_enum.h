// Copyright (c) 2013 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is a part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY, AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s) : Dmitry Anisimov, David Bommes, Kai Hormann, and Pierre Alliez.

/*!
  \file barycentric_enum.h
*/

#ifndef CGAL_BARYCENTRIC_ENUM_H
#define CGAL_BARYCENTRIC_ENUM_H

// CGAL namespace.
namespace CGAL {

/*!
 * \ingroup PkgBarycentric_coordinates_2
 * Barycentric_coordinates namespace contains implementations of all the generalized barycentric
 * coordinates: 2D, 3D, related enumerations, and so on.
 * It has an alies (short name) called `BC`.
 */

namespace Barycentric_coordinates {

// Possible locations of the query point provided by the user.
// UNSPECIFIED_LOCATION - Location is not known apriori, and it is defined automatically by the algorithm.
// AT_VERTEX - Query point is located at one of the polygon's vertices.
// ON_BOUNDARY - Query point is located on the boundary of the polygon.
// ON_BOUNDED_SIDE - Query point is located in the interior of the polygon excluding boundary.
// ON_UNBOUNDED_SIDE - Query point is located in the exterior of the polygon.

/// \name Locations of a query point
/// @{

/// Query_point_location is enumeration with possible locations of a
/// query point provided by the user.
enum Query_point_location
{
    /// Location is not known apriori, and it is defined automatically by the algorithm.
    UNSPECIFIED_LOCATION,

    /// Query point is located at one of the polygon's vertices.
    AT_VERTEX,

    /// Query point is located on the boundary of the polygon.
    ON_BOUNDARY,

    /// Query point is located inside the polygon excluding boundary.
    ON_BOUNDED_SIDE,

    /// Query point is located outside the polygon excluding boundary.
    ON_UNBOUNDED_SIDE
};

/// @}

// We can have two different algorithms to compute coordinates.
// PRECISE - Default slow algorithm, which is as precise as possible.
// FAST - Fast algorithm, which is less precise but much faster.

/// \name Types of an algorithm
/// @{

/// Type_of_algorithm is enumeration with possible algorithms to compute coordinates.
enum Type_of_algorithm
{
    /// Default slow algorithm, which is as precise as possible.
    PRECISE,

    /// Fast algorithm, which is less precise but much faster.
    FAST
};

/// @}

// Possible types of the input polygon.
// CONCAVE - Concave polygon = Non-convex polygon.
// WEAKLY_CONVEX - This is a convex polygon with collinear vertices.
// STRICTLY_CONVEX - This is a convex polygon without collinear vertices.

/// \name Types of a polygon
/// @{

/// Type_of_polygon is enumeration with possible types of the input polygon.
/// It is used internally to precondition coordinates.
enum Type_of_polygon
{
    /// Concave polygon = non-convex polygon.
    CONCAVE,

    /// This is a convex polygon with collinear vertices.
    WEAKLY_CONVEX,

    /// This is a convex polygon without collinear vertices.
    STRICTLY_CONVEX
};

} // namespace Barycentric_coordinates

/// @}

// Create a short alias for the Barycentric_coordinates namespace.

/// \name Namespace alias
/// @{

/// A short name (alias) of the namespace `CGAL::Barycentric_coordinates`.
namespace BC = Barycentric_coordinates;

/// @}

} // namespace CGAL

#endif // CGAL_BARYCENTRIC_ENUM_H