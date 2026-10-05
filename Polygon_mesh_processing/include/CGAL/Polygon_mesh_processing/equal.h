// Copyright (c) 2026 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Valque Léo

#ifndef CGAL_POLYGON_MESH_PROCESSING_EQUAL_H
#define CGAL_POLYGON_MESH_PROCESSING_EQUAL_H

#include <boost/graph/graph_traits.hpp>
#include <CGAL/Polygon_mesh_processing/polygon_mesh_to_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>

#include <CGAL/boost/graph/named_params_helper.h>

#include <algorithm>
#include <vector>

namespace PMP = CGAL::Polygon_mesh_processing;

template <typename Point>
struct Point_less
{
  bool operator()(const Point& p, const Point& q) const
  {
    if(p < q) return true;
    if(q < p) return false;
    return false;
  }
};

namespace CGAL {

namespace Polygon_mesh_processing {

/**
 * Reorders the points and polygons in a canonical way.
 */
template <typename PointRange, typename PolygonRange>
bool
canonical_soup(PointRange& points,
               PolygonRange& polygons)
{
  using Point = typename PointRange::value_type;
  using Polygon = typename PolygonRange::value_type;
  using Kernel = typename Kernel_traits<Point>::Kernel;
  typename Kernel::Less_xyz_3 less_xyz_3 = Kernel().less_xyz_3_object();

  std::vector<std::size_t> new_to_old_index(points.size());
  std::iota(new_to_old_index.begin(), new_to_old_index.end(), std::size_t(0));
  std::sort(new_to_old_index.begin(), new_to_old_index.end(),
            [&](std::size_t i, std::size_t j)
            { return less_xyz_3(points[i], points[j]); });

  std::vector<std::size_t> old_to_new_index(points.size());
  for(std::size_t new_index = 0; new_index < points.size(); ++new_index)
    old_to_new_index[new_to_old_index[new_index]] = new_index;

  std::vector<Point> reordered_points(points.size());
  for(std::size_t new_index = 0; new_index < points.size(); ++new_index)
    reordered_points[new_index] = points[new_to_old_index[new_index]];
  points = std::move(reordered_points);

  for(Polygon& polygon : polygons)
  {
    Polygon face;
    face.reserve(polygon.size());
    for(std::size_t i : polygon)
      face.push_back(old_to_new_index[i]);

    auto first = std::min_element(face.begin(), face.end());
    std::rotate(face.begin(), first, face.end());
    polygon = std::move(face);
  }

  auto lexicographical_less = [](const Polygon& p1, const Polygon& p2)
  {
    return std::lexicographical_compare(p1.begin(), p1.end(),
                                        p2.begin(), p2.end());
  };
  std::sort(polygons.begin(), polygons.end(), lexicographical_less);
  return true;
}

/**
 * Checks whether two polygon soups are equal, ignoring the order of their
 * points and polygons. The input ranges are reordered in place.
 */
template <typename PointRange1, typename PolygonRange1, typename PointRange2, typename PolygonRange2>
bool are_equal_polygon_soups(
    PointRange1& points1,
    PolygonRange1& polygons1,
    PointRange2& points2,
    PolygonRange2& polygons2)
{
  if(polygons1.size() != polygons2.size() || points1.size() != points2.size())
    return false;
    
  // The number of distinct geometric points must also agree.
  canonical_soup(points1, polygons1);
  canonical_soup(points2, polygons2);

  return std::equal(polygons1.begin(), polygons1.end(), polygons2.begin(), polygons2.end()) &&
         std::equal(points1.begin(), points1.end(), points2.begin(), points2.end());
}

/**
 * Checks whether two polygon soups are equal, ignoring the order of their
 * points and polygons. The input ranges are not modified.
 */
template <typename PointRange1, typename PolygonRange1, typename PointRange2, typename PolygonRange2>
bool are_equal_polygon_soups(
    const PointRange1& points1,
    const PolygonRange1& polygons1,
    const PointRange2& points2,
    const PolygonRange2& polygons2)
{
  PointRange1 points1_copy(points1);
  PolygonRange1 polygons1_copy(polygons1);
  PointRange2 points2_copy(points2);
  PolygonRange2 polygons2_copy(polygons2);

  return are_equal_polygon_soups(points1_copy, polygons1_copy,
                                 points2_copy, polygons2_copy);
}

/**
 * Checks whether two polygon meshes are equal, ignoring the order of their
 * vertices and faces. The meshes are converted to polygon soups for the
 * comparison.
 */
template <typename Mesh1, typename Mesh2,
          typename NamedParameters1 = parameters::Default_named_parameters,
          typename NamedParameters2 = parameters::Default_named_parameters>
bool are_equal_meshes(const Mesh1& mesh1,
                      const Mesh2& mesh2,
                      const NamedParameters1& np1 = parameters::default_values(),
                      const NamedParameters2& np2 = parameters::default_values())
{
  using VPM_helper_1 = GetVertexPointMap<Mesh1, NamedParameters1>;
  using VPM_helper_2 = GetVertexPointMap<Mesh2, NamedParameters2>;
  using VPM1 = typename VPM_helper_1::const_type;
  using VPM2 = typename VPM_helper_2::const_type;
  VPM1 vpm1 = VPM_helper_1::get_const_map(np1, mesh1);
  VPM2 vpm2 = VPM_helper_2::get_const_map(np2, mesh2);
  using Point1 = typename boost::property_traits<VPM1>::value_type;
  using Point2 = typename boost::property_traits<VPM2>::value_type;

  std::vector<Point1> points1;
  std::vector<std::vector<std::size_t> > polygons1;
  polygon_mesh_to_polygon_soup(mesh1, points1, polygons1, parameters::vertex_point_map(vpm1));

  std::vector<Point2> points2;
  std::vector<std::vector<std::size_t> > polygons2;
  polygon_mesh_to_polygon_soup(mesh2, points2, polygons2, parameters::vertex_point_map(vpm2));

  return are_equal_polygon_soups(points1, polygons1,
                                points2, polygons2);
}


} // namespace Polygon_mesh_processing
} // namespace CGAL
#endif //CGAL_POLYGON_MESH_PROCESSING_EQUAL_H
