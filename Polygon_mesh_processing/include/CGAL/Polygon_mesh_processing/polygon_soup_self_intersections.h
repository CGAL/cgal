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
#include <CGAL/tags.h>
#include <CGAL/type_traits.h>

#include <boost/iterator/transform_iterator.hpp>

#include <algorithm>
#include <iterator>
#include <vector>

namespace CGAL::Polygon_mesh_processing {

template <class ConcurrencyTag = Sequential_tag,
          class PointRange,
          class PolygonRange,
          class CGAL_NP_TEMPLATE_PARAMETERS>
bool does_polygon_soup_self_intersect(const PointRange& points,
                                      PolygonRange polygons,
                                      const CGAL_NP_CLASS& np = parameters::default_values())
{
  using std::cbegin;
  using std::cend;
  using std::size;

  using Default_pm = typename GetPointMap<PointRange, CGAL_NP_CLASS>::const_type;
  auto pm = parameters::choose_parameter<Default_pm>(parameters::get_parameter(np, internal_np::point_map));
  using Point_3 = CGAL::cpp20::remove_cvref_t<decltype(get(pm, *cbegin(points)))>;

  //let's be a bit more permissive and allow duplicated vertices
  auto to_point =[pm](const auto& v) { return get(pm, v); };
  auto first = boost::make_transform_iterator(cbegin(points), to_point);
  auto last = boost::make_transform_iterator(cend(points), to_point);

  std::vector<Point_3> unique_points(first, last);

  merge_duplicate_points_in_polygon_soup(unique_points, polygons);

  bool is_pure_triangles = std::all_of(CGAL_MAYBE_EXEC_POLICY(ConcurrencyTag)
                                       cbegin(polygons), cend(polygons),
                                       [](const auto&f) { return size(f) == 3; });
  if(is_pure_triangles)
  {
    // if the polygon soup is pure triangles,
    // we can use the triangle soup self-intersection function
    return does_triangle_soup_self_intersect<ConcurrencyTag>(unique_points, polygons);
  }

  // otherwise, we need to triangulate the polygons beforehand
  using Id = CGAL::cpp20::remove_cvref_t<decltype(*cbegin(*cbegin(polygons)))>;
  using New_polygon_type = std::vector<Id>;

  auto to_std_vector = [](const auto& poly) { return New_polygon_type{cbegin(poly), cend(poly)}; };

  std::vector<New_polygon_type> triangles(boost::make_transform_iterator(cbegin(polygons), to_std_vector),
                                          boost::make_transform_iterator(cend(polygons), to_std_vector));
  bool OK = triangulate_polygons(points, triangles, np);

  if (!OK) return true;

  return does_triangle_soup_self_intersect<ConcurrencyTag>(unique_points, triangles, np);
}

} // namespace CGAL::Polygon_mesh_processing


#endif // DOXYGEN_RUNNING

#endif // CGAL_POLYGON_MESH_PROCESSING_POLYGON_SOUP_SELF_INTERSECTIONS_H
