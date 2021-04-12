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

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_POLYLINE_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_POLYLINE_H

#include <CGAL/license/Shape_detection.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing_on_polyline/One_ring_neighbor_query.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_polyline/Least_squares_line_fit_region.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_polyline/Least_squares_line_fit_sorting.h>

namespace CGAL {
namespace Shape_detection {
namespace internal {

template<
typename InputRange,
typename OutputIterator,
typename NamedParameters>
OutputIterator region_growing_polylines(
  const InputRange& polyline, OutputIterator regions, const NamedParameters& np) {

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

template<
typename InputRange,
typename OutputIterator>
OutputIterator region_growing_polylines(
  const InputRange& polyline, OutputIterator regions) {

  return region_growing_polylines(
    polyline, regions, CGAL::parameters::all_default());
}

} // namespace internal
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_POLYLINE_H
