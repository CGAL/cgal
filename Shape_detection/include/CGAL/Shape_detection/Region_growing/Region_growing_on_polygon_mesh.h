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

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_POLYGON_MESH_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_POLYGON_MESH_H

#include <CGAL/license/Shape_detection.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing_on_polygon_mesh/Polyline_graph.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_polygon_mesh/One_ring_neighbor_query.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_polygon_mesh/Least_squares_plane_fit_region.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_polygon_mesh/Least_squares_plane_fit_sorting.h>

namespace CGAL {
namespace Shape_detection {
namespace internal {

template<
typename GeomTraits,
typename OutputIterator,
typename NamedParameters>
OutputIterator region_growing_planes(
  const CGAL::Polyhedron_3<GeomTraits, CGAL::Polyhedron_items_3, CGAL::HalfedgeDS_vector>& polyhedron,
  OutputIterator regions, const NamedParameters& np) {

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
typename GeomTraits,
typename OutputIterator>
OutputIterator region_growing_planes(
  const CGAL::Polyhedron_3<GeomTraits, CGAL::Polyhedron_items_3, CGAL::HalfedgeDS_vector>& polyhedron,
  OutputIterator regions) {

  return region_growing_planes(
    polyhedron, regions, CGAL::parameters::all_default());
}

template<
typename PointType,
typename OutputIterator,
typename NamedParameters>
OutputIterator region_growing_planes(
  const CGAL::Surface_mesh<PointType>& surface_mesh, OutputIterator regions, const NamedParameters& np) {

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
typename PointType,
typename OutputIterator>
OutputIterator region_growing_planes(
  const CGAL::Surface_mesh<PointType>& surface_mesh, OutputIterator regions) {

  return region_growing_planes(
    surface_mesh, regions, CGAL::parameters::all_default());
}

} // namespace internal
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_POLYGON_MESH_H
