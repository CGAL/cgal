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

#ifndef CGAL_ISOSURFACING_3_INTERNAL_IMPLICIT_CARTESIAN_GRID_GEOMETRY_3_H
#define CGAL_ISOSURFACING_3_INTERNAL_IMPLICIT_CARTESIAN_GRID_GEOMETRY_3_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/internal/Grid_topology_3.h>

namespace CGAL {
namespace Isosurfacing {
namespace internal {

// Describes the geometry of a regular Cartesian grid.
// Positions are not stored but calculated on-the-fly from an offset and grid spacing.
template <class GeomTraits>
class Implicit_Cartesian_grid_geometry_3
{
public:
  using Geom_traits = GeomTraits;
  using Point_3 = typename Geom_traits::Point_3;
  using Vector_3 = typename Geom_traits::Vector_3;

  using Vertex_descriptor = typename Grid_topology_3::Vertex_descriptor;

private:
  const Vector_3 m_offset;
  const Vector_3 m_spacing;
  Geom_traits m_gt;

public:
  Implicit_Cartesian_grid_geometry_3(const Vector_3& offset,
                                     const Vector_3& spacing,
                                     const Geom_traits& gt = Geom_traits())
    : m_offset{offset},
      m_spacing{spacing},
      m_gt{gt}
  { }

  // gets the position of vertex v
  Point_3 operator()(const Vertex_descriptor& v) const
  {
    typename Geom_traits::Compute_x_3 x_coord = m_gt.compute_x_3_object();
    typename Geom_traits::Compute_y_3 y_coord = m_gt.compute_y_3_object();
    typename Geom_traits::Compute_z_3 z_coord = m_gt.compute_z_3_object();
    typename Geom_traits::Construct_point_3 point = m_gt.construct_point_3_object();

    return point(v[0] * x_coord(m_spacing) + x_coord(m_offset),
                 v[1] * y_coord(m_spacing) + y_coord(m_offset),
                 v[2] * z_coord(m_spacing) + z_coord(m_offset));
  }
};

} // namespace internal
} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_INTERNAL_IMPLICIT_CARTESIAN_GRID_GEOMETRY_3_H
