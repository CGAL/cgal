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
// Author(s)     : Dmitry Anisimov
//

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_POLYLINE_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_POLYLINE_H

#include <CGAL/license/Shape_detection.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing_on_polyline/One_ring_neighbor_query.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_polyline/Least_squares_line_fit_region.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_polyline/Least_squares_line_fit_sorting.h>

namespace CGAL {
namespace Shape_detection {

template<
typename Point_d,
typename OutputIterator,
typename NamedParameters>
OutputIterator region_growing_polylines(
  const std::vector<Point_d>& /* polyline */, OutputIterator regions, const NamedParameters& /* np */) {

  // using GeomTraits  = typename Kernel_traits<Point_d>::Kernel;
  // using Input_range = std::vector<Point_d>;
  // using Normal_map  = CGAL::Identity_property_map<Point_d>;

  // using Neighbor_query = Polyline::One_ring_neighbor_query<GeomTraits, Input_range, Point_map>;
  // using Region_type    = Polyline::Least_squares_line_fit_region<GeomTraits, Input_range, Point_map>;
  // using Sorting        = Polyline::Least_squares_line_fit_sorting<GeomTraits, Input_range, Neighbor_query, Point_map>;
  // using Region_growing = Region_growing<Input_range, Neighbor_query, Region_type, typename Sorting::Seed_map>;

  // Neighbor_query neighbor_query(polyline, np, Point_map());
  // Region_type region_type(polyline, np, Point_map());
  // Sorting sorting(polyline, neighbor_query, Point_map());
  // sorting.sort();

  // Region_growing region_growing(
  //   polyline, neighbor_query, region_type, sorting.seed_map());
  // region_growing.detect(regions);
  return regions;
}

template<
typename Point_d,
typename OutputIterator>
OutputIterator region_growing_polylines(
  const std::vector<Point_d>& polyline, OutputIterator regions) {

  return region_growing_polylines(
    polyline, regions, CGAL::parameters::all_default());
}

} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_POLYLINE_H
