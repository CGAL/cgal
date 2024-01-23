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

#ifndef CGAL_ISOSURFACING_3_INTERNAL_EXPLICIT_CARTESIAN_GRID_FUNCTION_H
#define CGAL_ISOSURFACING_3_INTERNAL_EXPLICIT_CARTESIAN_GRID_FUNCTION_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/internal/Grid_topology_3.h>

namespace CGAL {
namespace Isosurfacing {
namespace internal {

template <typename Grid_>
class Explicit_Cartesian_grid_function
{
public:
  using Grid = Grid_;
  using Geom_traits = typename Grid::Geom_traits;
  using FT = typename Geom_traits::FT;

  using Vertex_descriptor = typename Grid_topology_3::Vertex_descriptor;

public:
  Explicit_Cartesian_grid_function(const Grid& grid)
    : m_grid{grid}
  { }

  // gets the value at vertex `v`
  FT operator()(const Vertex_descriptor& v) const
  {
    return m_grid.value(v[0], v[1], v[2]);
  }

private:
  const Grid& m_grid;
};

} // namespace internal
} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_INTERNAL_EXPLICIT_CARTESIAN_GRID_FUNCTION_H
