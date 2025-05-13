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
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>

#include <boost/iterator/transform_iterator.hpp>

namespace CGAL::Polygon_mesh_processing {

template <class ConcurrencyTag = Sequential_tag,
          class PointRange,
          class PolygonRange,
          class CGAL_NP_TEMPLATE_PARAMETERS>
bool does_polygon_soup_self_intersect(const PointRange& points,
                                      PolygonRange polygons,
                                      const CGAL_NP_CLASS& np = parameters::default_values())
{
  //let's be a bit more permissive and allow duplicated vertices
  auto pm = parameters::choose_parameter<typename GetPointMap<PointRange, CGAL_NP_CLASS>::const_type>(parameters::get_parameter(np, internal_np::point_map));
  using Point_3 = typename boost::property_traits<decltype(pm)>::value_type;
  using PointRange_const_iterator = typename PointRange::const_iterator;
  using PointRange_value_type = typename std::iterator_traits<PointRange_const_iterator>::value_type;
  auto to_point =[pm](const PointRange_value_type& v) { return get(pm, v); };

  std::vector<Point_3> unique_points(boost::make_transform_iterator(points.begin(), to_point),
                                     boost::make_transform_iterator(points.end(), to_point));

  merge_duplicate_points_in_polygon_soup(unique_points, polygons);

  bool is_pure_triangles = std::all_of(polygons.begin(), polygons.end(),
                                       [](const auto&f) { return f.size() == 3; });
  if(is_pure_triangles)
  {
    // if the polygon soup is pure triangles,
    // we can use the triangle soup self-intersection function
    return does_triangle_soup_self_intersect<ConcurrencyTag>(unique_points, polygons);
  }

  // otherwise, we need to triangulate the polygons beforehand
  using Polygon = CGAL::cpp20::remove_cvref_t<decltype(*polygons.begin())>;
  using Id = typename std::iterator_traits<typename Polygon::const_iterator>::value_type;
  auto to_std_vector = [](const Polygon& poly)
  {
    return std::vector<Id>(poly.begin(), poly.end());
  };

  std::vector<std::vector<Id>> triangles(boost::make_transform_iterator(polygons.begin(), to_std_vector),
                                         boost::make_transform_iterator(polygons.end(), to_std_vector));
  bool OK = triangulate_polygons(points, triangles, np);

  if (!OK) return false;

  return does_triangle_soup_self_intersect<ConcurrencyTag>(unique_points, triangles, np);
}

} // namespace CGAL::Polygon_mesh_processing


#endif // DOXYGEN_RUNNING

#endif // CGAL_POLYGON_MESH_PROCESSING_POLYGON_SOUP_SELF_INTERSECTIONS_H
