// Copyright (c) 2023  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_CDT_3_DEBUG_IO_H
#define CGAL_CDT_3_DEBUG_IO_H

#include <ostream>
#include <CGAL/boost/graph/IO/OFF.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

namespace CGAL {


  template <typename Tr, typename Facets>
  auto export_facets_to_surface_mesh(const Tr& tr, Facets&& facets_range) {
    using Point_3 = typename Tr::Geom_traits::Point_3;
    const auto size = std::distance(facets_range.begin(), facets_range.end());
    std::vector<Point_3> points;
    points.reserve(size * 3);
    std::vector<std::array<std::size_t, 3>> facets;
    facets.reserve(size);

    for(std::size_t i = 0; const auto& [cell, facet_index] : facets_range) {
      const auto v0 = cell->vertex(Tr::vertex_triple_index(facet_index, 0));
      const auto v1 = cell->vertex(Tr::vertex_triple_index(facet_index, 1));
      const auto v2 = cell->vertex(Tr::vertex_triple_index(facet_index, 2));
      points.push_back(tr.point(v0));
      points.push_back(tr.point(v1));
      points.push_back(tr.point(v2));
      facets.push_back({i, i+1, i+2});
      i += 3;
    }
    CGAL::Polygon_mesh_processing::merge_duplicate_points_in_polygon_soup(points, facets);
    CGAL::Polygon_mesh_processing::orient_polygon_soup(points, facets);
    Surface_mesh<Point_3> mesh;
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, facets, mesh);
    return mesh;
  }

  template <typename Tr, typename Facets>
  void write_facets(std::ostream& out, const Tr& tr, Facets&& facets_range) {
    const auto mesh = export_facets_to_surface_mesh(tr, std::forward<Facets>(facets_range));
    CGAL::IO::write_OFF(out, mesh);
  }

} // end namespace CGAL

#endif // CGAL_CDT_3_DEBUG_IO_H
