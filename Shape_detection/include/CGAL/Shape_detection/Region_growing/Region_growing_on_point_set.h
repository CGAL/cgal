// Copyright (c) 2018 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Florent Lafarge, Simon Giraudot, Thien Hoang, Dmitry Anisimov
//

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_POINT_SET_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_POINT_SET_H

#include <CGAL/license/Shape_detection.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing_on_point_set/K_neighbor_query.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_point_set/Sphere_neighbor_query.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing_on_point_set/Least_squares_line_fit_region.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_point_set/Least_squares_plane_fit_region.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing_on_point_set/Least_squares_line_fit_sorting.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_point_set/Least_squares_plane_fit_sorting.h>

#include <CGAL/Point_set_3.h>

namespace CGAL {
namespace Shape_detection {
namespace internal {

template<
typename InputRange,
typename OutputIterator,
typename NamedParameters>
OutputIterator region_growing_lines(
  const InputRange& points_with_normals, OutputIterator regions, const NamedParameters& np) {

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
typename OutputIterator>
OutputIterator region_growing_lines(
  const InputRange& points_with_normals, OutputIterator regions) {

  return region_growing_lines(
    points_with_normals, regions, CGAL::parameters::all_default());
}

template<
typename InputRange,
typename OutputIterator,
typename NamedParameters>
OutputIterator region_growing_planes(
  const InputRange& points_with_normals, OutputIterator regions, const NamedParameters& np) {

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
typename InputRange,
typename OutputIterator>
OutputIterator region_growing_planes(
  const InputRange& points_with_normals, OutputIterator regions) {

  return region_growing_planes(
    points_with_normals, regions, CGAL::parameters::all_default());
}

template<
typename PointType,
typename VectorType,
typename OutputIterator,
typename NamedParameters>
OutputIterator region_growing_planes(
  const CGAL::Point_set_3<PointType, VectorType>& point_set, OutputIterator regions, const NamedParameters& np) {

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
typename OutputIterator,
typename NamedParameters>
OutputIterator region_growing_planes(
  const CGAL::Point_set_3<PointType, VectorType>& point_set, OutputIterator regions) {

  return region_growing_planes(
    point_set, regions, CGAL::parameters::all_default());
}

} // namespace internal
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_POINT_SET_H
