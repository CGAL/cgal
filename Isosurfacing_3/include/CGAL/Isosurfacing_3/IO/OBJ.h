// Copyright (c) 2022-2024 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_ISOSURFACING_3_IMAGE_3_H
#define CGAL_ISOSURFACING_3_IMAGE_3_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/Cartesian_grid_3.h>

#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Named_function_parameters.h>

#include <string>
#include <fstream>

namespace CGAL {
namespace Isosurfacing {

template <typename GeomTraits, typename MemoryPolicy>
class Cartesian_grid_3;

namespace IO {

template <typename GeomTraits, typename MemoryPolicy,
          typename NamedParameters = parameters::Default_named_parameters>
bool write_OBJ(std::ostream& out,
               const Cartesian_grid_3<GeomTraits, MemoryPolicy>& grid,
               const NamedParameters& np = parameters::default_values())
{
  using Point_3 = typename GeomTraits::Point_3;

  auto x_coord = grid.geom_traits().compute_x_3_object();
  auto y_coord = grid.geom_traits().compute_y_3_object();
  auto z_coord = grid.geom_traits().compute_z_3_object();
  auto vertex = grid.geom_traits().construct_vertex_3_object();

  set_ascii_mode(out); // obj is ASCII only

  set_stream_precision_from_NP(out, np);

  if(out.fail())
    return false;

  // write vertices
  for(std::size_t x=0; x<grid.xdim(); ++x) {
    for(std::size_t y=0; y<grid.ydim(); ++y) {
      for(std::size_t z=0; z<grid.zdim(); ++z)
      {
        const Point_3& p = vertex(grid.span(), 0);
        const double x_coord_d = CGAL::to_double(x_coord(p) + x * grid.spacing()[0]);
        const double y_coord_d = CGAL::to_double(y_coord(p) + y * grid.spacing()[1]);
        const double z_coord_d = CGAL::to_double(z_coord(p) + z * grid.spacing()[2]);
        out << "v " << x_coord_d << " " << y_coord_d << " " << z_coord_d << std::endl;
      }
    }
  }

  // write faces
  for(std::size_t x=0; x<grid.xdim()-1; ++x) {
    for(std::size_t y=0; y<grid.ydim()-1; ++y) {
      for(std::size_t z=0; z<grid.zdim()-1; ++z)
      {
        const std::size_t v0 = (z * grid.ydim() + y) * grid.xdim() + x;
        const std::size_t v1 = (z * grid.ydim() + y + 1) * grid.xdim() + x;
        const std::size_t v2 = (z * grid.ydim() + y + 1) * grid.xdim() + x + 1;
        const std::size_t v3 = (z * grid.ydim() + y) * grid.xdim() + x + 1;
        out << "f " << v0+1 << " " << v1+1 << " " << v2+1 << " " << v3+1 << std::endl;
      }
    }
  }

  return out.good();
}

template <typename GeomTraits, typename MemoryPolicy,
          typename NamedParameters = parameters::Default_named_parameters>
bool write_OBJ(const std::string& fname,
               const Cartesian_grid_3<GeomTraits, MemoryPolicy>& grid,
               const NamedParameters& np = parameters::default_values())
{
  std::ofstream os(fname);
  CGAL::IO::set_mode(os, CGAL::IO::ASCII);

  return write_OBJ(os, grid, np);
}

} // namespace IO
} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_IMAGE_3_H
