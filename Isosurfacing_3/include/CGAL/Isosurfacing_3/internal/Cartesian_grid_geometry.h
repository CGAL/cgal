// Copyright (c) 2022-2023 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Julian Stahl

#ifndef CGAL_ISOSURFACING_3_INTERNAL_CARTESIAN_GRID_GEOMETRY_H
#define CGAL_ISOSURFACING_3_INTERNAL_CARTESIAN_GRID_GEOMETRY_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/internal/Grid_topology.h>

namespace CGAL {
namespace Isosurfacing {
namespace internal {

// describes the geometry of a regular Cartesian grid.
// Positions are not stored but calculated from an offset and grid spacing.
template <typename GeomTraits>
class Cartesian_grid_geometry
{
public:
  using Geom_traits = GeomTraits;
  using Point = typename Geom_traits::Point_3;
  using Vector = typename Geom_traits::Vector_3;

  using Vertex_descriptor = typename Grid_topology::Vertex_descriptor;

public:
  // creates a regular grid geometry where `offset` is the position of the vertex with index `(0, 0, 0)`
  // and `spacing` the distance between two connected vertices in each dimension.
  Cartesian_grid_geometry(const Vector& offset,
                          const Vector& spacing)
    : offset(offset),
      spacing(spacing)
  { }

  // gets the position of vertex `v`
  Point operator()(const Vertex_descriptor& v) const
  {
    return Point(v[0] * spacing[0],
                 v[1] * spacing[1],
                 v[2] * spacing[2]) + offset;
  }

private:
  Vector offset;
  Vector spacing;
};

} // namespace internal
} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_INTERNAL_CARTESIAN_GRID_GEOMETRY_H
