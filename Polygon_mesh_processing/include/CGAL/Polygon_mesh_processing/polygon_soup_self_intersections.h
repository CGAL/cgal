// Copyright (c) 2025 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_POLYGON_MESH_PROCESSING_POLYGON_SOUP_SELF_INTERSECTIONS_H
#define CGAL_POLYGON_MESH_PROCESSING_POLYGON_SOUP_SELF_INTERSECTIONS_H

#include <CGAL/license/Polygon_mesh_processing/predicate.h>

#ifndef DOXYGEN_RUNNING

#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

namespace CGAL::Polygon_mesh_processing {

template <class ConcurrencyTag = Sequential_tag,
          class PointRange,
          class PolygonRange,
          class CGAL_NP_TEMPLATE_PARAMETERS>
bool does_polygon_soup_self_intersect(const PointRange& points,
                                      const PolygonRange& polygons,
                                      const CGAL_NP_CLASS& np = parameters::default_values())
{
  bool is_pure_triangles = std::all_of(polygons.begin(), polygons.end(),
                                       [](const auto&f) { return f.size() == 3; });
  if(is_pure_triangles)
  {
    // if the polygon soup is pure triangles,
    // we can use the triangle soup self-intersection function
    return does_triangle_soup_self_intersect<ConcurrencyTag>(points, polygons, np);
  }

  // otherwise, we need to triangulate the polygons beforehand
  using Polygon = CGAL::cpp20::remove_cvref_t<decltype(*polygons.begin())>;
  std::vector<Polygon> triangles(polygons.begin(), polygons.end());
  triangulate_polygons(points, triangles, np);

  return does_triangle_soup_self_intersect<ConcurrencyTag>(points, triangles, np);
}

} // namespace CGAL::Polygon_mesh_processing


#endif // DOXYGEN_RUNNING

#endif // CGAL_POLYGON_MESH_PROCESSING_POLYGON_SOUP_SELF_INTERSECTIONS_H
