// Copyright (c) 2020 GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Dmitry Anisimov
//

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_INTERNAL_FREE_FUNCTIONS_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_INTERNAL_FREE_FUNCTIONS_H

#include <CGAL/license/Shape_detection.h>

#include <CGAL/Point_set_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_point_set.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_polygon_mesh.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_polyline.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_segment_set.h>

namespace CGAL {
namespace Shape_detection {
namespace internal {

template<
typename InputRange,
typename OutputIterator,
typename CGAL_NP_TEMPLATE_PARAMETERS>
OutputIterator region_growing_lines(
  const InputRange& points_with_normals, OutputIterator regions, const CGAL_NP_CLASS& np = parameters::default_values()) {

  using Input_range = InputRange;
  using Point_with_normal = typename Input_range::value_type;
  using Point_type = typename Point_with_normal::first_type;
  using Traits = typename Kernel_traits<Point_type>::Kernel;
  using Point_map = CGAL::First_of_pair_property_map<Point_with_normal>;
  using Normal_map = CGAL::Second_of_pair_property_map<Point_with_normal>;

  using Neighbor_query = Point_set::K_neighbor_query<Traits, Input_range, Point_map>;
  using Region_type = Point_set::Least_squares_line_fit_region<Traits, Input_range, Point_map, Normal_map>;
  using Sorting = Point_set::Least_squares_line_fit_sorting<Traits, Input_range, Neighbor_query, Point_map>;
  using Region_growing = Region_growing<Input_range, Neighbor_query, Region_type, typename Sorting::Seed_map>;

  Neighbor_query neighbor_query(points_with_normals, np);
  Region_type region_type(points_with_normals, np);
  Sorting sorting(points_with_normals, neighbor_query, np);
  sorting.sort();

  Region_growing region_growing(
    points_with_normals, neighbor_query, region_type, sorting.seed_map());
  region_growing.detect(regions);
  return regions;
}

template<
typename InputRange,
typename OutputIterator,
typename CGAL_NP_TEMPLATE_PARAMETERS>
OutputIterator region_growing_planes(
  const InputRange& points_with_normals, OutputIterator regions, const CGAL_NP_CLASS& np = parameters::default_values()) {

  using Input_range = InputRange;
  using Point_with_normal = typename Input_range::value_type;
  using Point_type = typename Point_with_normal::first_type;
  using Traits = typename Kernel_traits<Point_type>::Kernel;
  using Point_map = CGAL::First_of_pair_property_map<Point_with_normal>;
  using Normal_map = CGAL::Second_of_pair_property_map<Point_with_normal>;

  using Neighbor_query = Point_set::K_neighbor_query<Traits, Input_range, Point_map>;
  using Region_type = Point_set::Least_squares_plane_fit_region<Traits, Input_range, Point_map, Normal_map>;
  using Sorting = Point_set::Least_squares_plane_fit_sorting<Traits, Input_range, Neighbor_query, Point_map>;
  using Region_growing = Region_growing<Input_range, Neighbor_query, Region_type, typename Sorting::Seed_map>;

  Neighbor_query neighbor_query(points_with_normals, np);
  Region_type region_type(points_with_normals, np);
  Sorting sorting(points_with_normals, neighbor_query, np);
  sorting.sort();

  Region_growing region_growing(
    points_with_normals, neighbor_query, region_type, sorting.seed_map());
  region_growing.detect(regions);
  return regions;
}

template<
typename PointType,
typename VectorType,
typename OutputIterator,
typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT>
OutputIterator region_growing_planes(
  const CGAL::Point_set_3<PointType, VectorType>& point_set, OutputIterator regions, const CGAL_NP_CLASS& np) {

  using Point_type = PointType;
  using Vector_type = VectorType;
  using Traits = typename Kernel_traits<Point_type>::Kernel;
  using Point_set_3 = CGAL::Point_set_3<Point_type, Vector_type>;
  using Point_map = typename Point_set_3::Point_map;
  using Normal_map = typename Point_set_3::Vector_map;

  using Neighbor_query = Point_set::K_neighbor_query<Traits, Point_set_3, Point_map>;
  using Region_type = Point_set::Least_squares_plane_fit_region<Traits, Point_set_3, Point_map, Normal_map>;
  using Sorting = Point_set::Least_squares_plane_fit_sorting<Traits, Point_set_3, Neighbor_query, Point_map>;
  using Region_growing = Region_growing<Point_set_3, Neighbor_query, Region_type, typename Sorting::Seed_map>;

  CGAL_precondition(point_set.has_normal_map());
  Neighbor_query neighbor_query(point_set, np);
  Region_type region_type(point_set, np);
  Sorting sorting(point_set, neighbor_query, np);
  sorting.sort();

  Region_growing region_growing(
    point_set, neighbor_query, region_type, sorting.seed_map());
  region_growing.detect(regions);
  return regions;
}

template<
typename PointType,
typename VectorType,
typename OutputIterator>
OutputIterator region_growing_planes(
  const CGAL::Point_set_3<PointType, VectorType>& point_set, OutputIterator regions) {

  return region_growing_planes(
    point_set, regions, CGAL::parameters::
    point_map(point_set.point_map()).
    normal_map(point_set.normal_map()));
}

template<
typename GeomTraits,
typename OutputIterator,
typename CGAL_NP_TEMPLATE_PARAMETERS>
OutputIterator region_growing_planes(
  const CGAL::Polyhedron_3<GeomTraits, CGAL::Polyhedron_items_3, CGAL::HalfedgeDS_vector>& polyhedron,
  OutputIterator regions, const CGAL_NP_CLASS& np = parameters::default_values()) {

  using Traits = GeomTraits;
  using Polyhedron = CGAL::Polyhedron_3<Traits, CGAL::Polyhedron_items_3, CGAL::HalfedgeDS_vector>;
  using Face_range = typename CGAL::Iterator_range<typename boost::graph_traits<Polyhedron>::face_iterator>;
  using Neighbor_query = Polygon_mesh::One_ring_neighbor_query<Polyhedron, Face_range>;
  using Region_type = Polygon_mesh::Least_squares_plane_fit_region<Traits, Polyhedron, Face_range>;
  using Sorting = Polygon_mesh::Least_squares_plane_fit_sorting<Traits, Polyhedron, Neighbor_query, Face_range>;
  using Region_growing = Region_growing<Face_range, Neighbor_query, Region_type, typename Sorting::Seed_map>;

  const Face_range face_range = faces(polyhedron);
  Neighbor_query neighbor_query(polyhedron);
  Region_type region_type(polyhedron, np);
  Sorting sorting(polyhedron, neighbor_query, np);
  sorting.sort();

  Region_growing region_growing(
    face_range, neighbor_query, region_type, sorting.seed_map());
  region_growing.detect(regions);
  return regions;
}

template<
typename PointType,
typename OutputIterator,
typename CGAL_NP_TEMPLATE_PARAMETERS>
OutputIterator region_growing_planes(
  const CGAL::Surface_mesh<PointType>& surface_mesh, OutputIterator regions, const CGAL_NP_CLASS& np = parameters::default_values()) {

  using Point_type = PointType;
  using Traits = typename Kernel_traits<Point_type>::Kernel;
  using Surface_mesh = CGAL::Surface_mesh<Point_type>;
  using Face_range = typename Surface_mesh::Face_range;
  using Neighbor_query = Polygon_mesh::One_ring_neighbor_query<Surface_mesh>;
  using Region_type = Polygon_mesh::Least_squares_plane_fit_region<Traits, Surface_mesh>;
  using Sorting = Polygon_mesh::Least_squares_plane_fit_sorting<Traits, Surface_mesh, Neighbor_query>;
  using Region_growing = Region_growing<Face_range, Neighbor_query, Region_type, typename Sorting::Seed_map>;

  const Face_range face_range = faces(surface_mesh);
  Neighbor_query neighbor_query(surface_mesh);
  Region_type region_type(surface_mesh, np);
  Sorting sorting(surface_mesh, neighbor_query, np);
  sorting.sort();

  Region_growing region_growing(
    face_range, neighbor_query, region_type, sorting.seed_map());
  region_growing.detect(regions);
  return regions;
}

template<
typename InputRange,
typename OutputIterator,
typename CGAL_NP_TEMPLATE_PARAMETERS>
OutputIterator region_growing_polylines(
  const InputRange& polyline, OutputIterator regions, const CGAL_NP_CLASS& np = parameters::default_values()) {

  using Input_range = InputRange;
  using Point_type = typename Input_range::value_type;
  using Traits = typename Kernel_traits<Point_type>::Kernel;
  using Point_map = CGAL::Identity_property_map<Point_type>;

  using Neighbor_query = Polyline::One_ring_neighbor_query<Traits, Input_range>;
  using Region_type = Polyline::Least_squares_line_fit_region<Traits, Input_range, Point_map>;
  using Sorting = Polyline::Least_squares_line_fit_sorting<Traits, Input_range, Neighbor_query, Point_map>;
  using Region_growing = Region_growing<Input_range, Neighbor_query, Region_type, typename Sorting::Seed_map>;

  Neighbor_query neighbor_query(polyline);
  Region_type region_type(polyline, np);
  Sorting sorting(polyline, neighbor_query, np);
  sorting.sort();

  Region_growing region_growing(
    polyline, neighbor_query, region_type, sorting.seed_map());
  region_growing.detect(regions);
  return regions;
}

} // namespace internal
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_INTERNAL_FREE_FUNCTIONS_H
