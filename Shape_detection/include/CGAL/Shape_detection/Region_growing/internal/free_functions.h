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

#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Point_set.h>
#include <CGAL/Shape_detection/Region_growing/Triangle_mesh.h>
#include <CGAL/Shape_detection/Region_growing/Polyline.h>
#include <CGAL/Shape_detection/Region_growing/Segment_set.h>

// TODO: this is not using generic programming for FaceGraph and PointSet
//       as a consequence, artificial dependencies on Surface_mesh, Polyhedron,
//       HDS and Point_set_3 are done.

namespace CGAL {
namespace Shape_detection {
namespace internal {

template<
typename InputRange,
typename OutputIterator,
typename CGAL_NP_TEMPLATE_PARAMETERS>
OutputIterator region_growing_lines(
  const InputRange& points, OutputIterator regions, const CGAL_NP_CLASS& np = parameters::default_values()) {

  // basic geometric types
  typedef Point_set_processing_3_np_helper<InputRange, CGAL_NP_CLASS> NP_helper;
  typedef typename NP_helper::Point_map PointMap;
  typedef typename NP_helper::Normal_map NormalMap;
  typedef typename NP_helper::Geom_traits Kernel;

  using Neighbor_query = Point_set::K_neighbor_query<Kernel, InputRange, PointMap>;
  using Region_type = Point_set::Least_squares_line_fit_region<Kernel, InputRange, PointMap, NormalMap>;
  using Sorting = Point_set::Least_squares_line_fit_sorting<Kernel, InputRange, Neighbor_query, PointMap>;
  using Region_growing = Region_growing<InputRange, Neighbor_query, Region_type, typename Sorting::Seed_map>;

  Neighbor_query neighbor_query(points, np);
  Region_type region_type(points, np);
  Sorting sorting(points, neighbor_query, np);
  sorting.sort();

  Region_growing region_growing(
    points, neighbor_query, region_type, sorting.seed_map());
  region_growing.detect(regions);
  return regions;
}

template<
typename InputRange,
typename OutputIterator,
typename CGAL_NP_TEMPLATE_PARAMETERS>
OutputIterator region_growing_planes(
  const InputRange& points, OutputIterator regions, const CGAL_NP_CLASS& np = parameters::default_values()) {

  // basic geometric types
  typedef Point_set_processing_3_np_helper<InputRange, CGAL_NP_CLASS> NP_helper;
  typedef typename NP_helper::Point_map PointMap;
  typedef typename NP_helper::Normal_map NormalMap;
  typedef typename NP_helper::Geom_traits Kernel;

  using Neighbor_query = Point_set::K_neighbor_query<Kernel, InputRange, PointMap>;
  using Region_type = Point_set::Least_squares_plane_fit_region<Kernel, InputRange, PointMap, NormalMap>;
  using Sorting = Point_set::Least_squares_plane_fit_sorting<Kernel, InputRange, Neighbor_query, PointMap>;
  using Region_growing = Region_growing<InputRange, Neighbor_query, Region_type, typename Sorting::Seed_map>;

  Neighbor_query neighbor_query(points, np);
  Region_type region_type(points, np);
  Sorting sorting(points, neighbor_query, np);
  sorting.sort();

  Region_growing region_growing(
    points, neighbor_query, region_type, sorting.seed_map());
  region_growing.detect(regions);
  return regions;
}

template<
typename TriangleMesh,
typename OutputIterator,
typename CGAL_NP_TEMPLATE_PARAMETERS>
OutputIterator region_growing_planes_triangle_mesh(
  const TriangleMesh& triangle_mesh, OutputIterator regions, const CGAL_NP_CLASS& np = parameters::default_values()) {

  using Kernel = typename Kernel_traits<typename TriangleMesh::Point>::Kernel;
  using Face_iterator = typename boost::graph_traits<TriangleMesh>::face_iterator;
  using Face_range = Iterator_range<Face_iterator>;

  using Neighbor_query = Triangle_mesh::One_ring_neighbor_query<TriangleMesh>;
  using Region_type = Triangle_mesh::Least_squares_plane_fit_region<Kernel, TriangleMesh, Face_range>;
  using Sorting = Triangle_mesh::Least_squares_plane_fit_sorting<Kernel, TriangleMesh, Neighbor_query, Face_range>;
  using Region_growing = Region_growing<Face_range, Neighbor_query, Region_type, typename Sorting::Seed_map>;

  Neighbor_query neighbor_query(triangle_mesh);
  Region_type region_type(triangle_mesh, np);
  Sorting sorting(triangle_mesh, neighbor_query, np);
  sorting.sort();

  Region_growing region_growing(
    faces(triangle_mesh), neighbor_query, region_type, sorting.seed_map());
  region_growing.detect(regions);
  return regions;
}

template<
typename InputRange,
typename OutputIterator,
typename CGAL_NP_TEMPLATE_PARAMETERS>
OutputIterator region_growing_polylines(
  const InputRange& polyline, OutputIterator regions, const CGAL_NP_CLASS& np = parameters::default_values()) {

  using Point_type = typename InputRange::value_type;
  using Kernel = typename Kernel_traits<Point_type>::Kernel;
  using Point_map = CGAL::Identity_property_map<Point_type>;

  using Neighbor_query = Polyline::One_ring_neighbor_query<Kernel, InputRange>;
  using Region_type = Polyline::Least_squares_line_fit_region<Kernel, InputRange, Point_map>;
  using Sorting = Polyline::Least_squares_line_fit_sorting<Kernel, InputRange, Neighbor_query, Point_map>;
  using Region_growing = Region_growing<InputRange, Neighbor_query, Region_type, typename Sorting::Seed_map>;

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
