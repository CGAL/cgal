// Copyright (c) 2019 GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_KSR_3_RECONSTRUCTION_H
#define CGAL_KSR_3_RECONSTRUCTION_H

// #include <CGAL/license/Kinetic_shape_reconstruction.h>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_point_set.h>

// Internal includes.
#include <CGAL/KSR/enum.h>
#include <CGAL/KSR/utils.h>
#include <CGAL/KSR/debug.h>
#include <CGAL/KSR/property_map.h>
#include <CGAL/KSR_3/Data_structure.h>
#include <CGAL/KSR_3/Visibility.h>
#include <CGAL/KSR_3/Graphcut.h>

namespace CGAL {
namespace KSR_3 {

template<
typename InputRange,
typename PointMap,
typename VectorMap,
typename SemanticMap,
typename GeomTraits>
class Reconstruction {

public:
  using Input_range  = InputRange;
  using Point_map    = PointMap;
  using Vector_map   = VectorMap;
  using Semantic_map = SemanticMap;
  using Kernel       = GeomTraits;

private:
  using FT      = typename Kernel::FT;
  using Point_2 = typename Kernel::Point_2;
  using Point_3 = typename Kernel::Point_3;
  using Plane_3 = typename Kernel::Plane_3;

  using Data_structure = KSR_3::Data_structure<Kernel>;
  using Point_map_3    = KSR::Item_property_map<Input_range, Point_map>;
  using Vector_map_3   = KSR::Item_property_map<Input_range, Vector_map>;
  using PFace          = typename Data_structure::PFace;

  using Semantic_label    = KSR::Semantic_label;
  using Planar_shape_type = KSR::Planar_shape_type;

  using Indices     = std::vector<std::size_t>;
  using Polygon_3   = std::vector<Point_3>;
  using Polygon_map = CGAL::Identity_property_map<Polygon_3>;

  using IK       = Exact_predicates_inexact_constructions_kernel;
  using IPoint_3  = typename IK::Point_3;
  using IPlane_3  = typename IK::Plane_3;

  using Converter = CGAL::Cartesian_converter<Kernel, IK>;
  using Delaunay  = CGAL::Delaunay_triangulation_2<Kernel>;

  // using Neighbor_query_3 = CGAL::Shape_detection::Point_set::
  //   Sphere_neighbor_query<Kernel, Indices, Point_map_3>;

  using Neighbor_query_3 = CGAL::Shape_detection::Point_set::
    K_neighbor_query<Kernel, Indices, Point_map_3>;
  using Planar_region    = CGAL::Shape_detection::Point_set::
    Least_squares_plane_fit_region<Kernel, Indices, Point_map_3, Vector_map_3>;
  using Planar_sorting   = CGAL::Shape_detection::Point_set::
    Least_squares_plane_fit_sorting<Kernel, Indices, Neighbor_query_3, Point_map_3>;
  using Region_growing   = CGAL::Shape_detection::
    Region_growing<Indices, Neighbor_query_3, Planar_region, typename Planar_sorting::Seed_map>;

  using Visibility_label = KSR::Visibility_label;
  using Visibility       = KSR_3::Visibility<Kernel, Point_map_3, Vector_map_3>;
  using Graphcut         = KSR_3::Graphcut<Kernel>;

public:

  Reconstruction(
    const Input_range& input_range,
    const Point_map& point_map,
    const Vector_map& normal_map,
    const Semantic_map& semantic_map,
    Data_structure& data,
    const bool verbose,
    const bool debug) :
  m_input_range(input_range),
  m_semantic_map(semantic_map),
  m_point_map_3(m_input_range, point_map),
  m_normal_map_3(m_input_range, normal_map),
  m_data(data),
  m_debug(debug),
  m_verbose(verbose),
  m_planar_shape_type(Planar_shape_type::CONVEX_HULL) {

    clear();
    collect_points(Semantic_label::GROUND           , m_ground_points);
    collect_points(Semantic_label::BUILDING_BOUNDARY, m_boundary_points);
    collect_points(Semantic_label::BUILDING_INTERIOR, m_interior_points);

    if (m_verbose) {
      std::cout << std::endl << "--- RECONSTRUCTION: " << std::endl;
      std::cout << "* num ground points: "   << m_ground_points.size()   << std::endl;
      std::cout << "* num boundary points: " << m_boundary_points.size() << std::endl;
      std::cout << "* num interior points: " << m_interior_points.size() << std::endl;
    }
  }

  template<typename NamedParameters>
  const bool detect_planar_shapes(
    const NamedParameters& np) {

    if (m_verbose) {
      std::cout << std::endl << "--- DETECTING PLANAR SHAPES: " << std::endl;
    }
    m_polygons.clear();
    create_ground_plane();
    CGAL_assertion(m_polygons.size() == 1);
    create_approximate_walls(np);
    create_approximate_roofs(np);
    if (m_debug) dump_polygons("detected-planar-shapes");
    return true;
  }

  template<typename NamedParameters>
  const bool regularize_planar_shapes(
    const NamedParameters& np) {

    return true;
    CGAL_assertion_msg(false, "TODO: REGULARIZE PLANAR SHAPES!");
    return false;
  }

  template<typename NamedParameters>
  const bool compute_model(
    const NamedParameters& np) {

    if (m_verbose) {
      std::cout << std::endl << "--- COMPUTING THE MODEL: " << std::endl;
    }

    std::map<PFace, Indices> pface_points;
    assign_points_to_pfaces(pface_points);
    const Visibility visibility(
      m_data, pface_points, m_point_map_3, m_normal_map_3);

    CGAL_assertion(m_data.volumes().size() > 0);
    visibility.compute(m_data.volumes());
    if (m_debug) dump_volumes("visibility");

    const FT beta = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::graphcut_beta), FT(1) / FT(2));

    Graphcut graphcut(m_data, beta);
    graphcut.compute(m_data.volumes());
    if (m_debug) dump_volumes("graphcut");

    extract_surface();
    CGAL_assertion_msg(false, "TODO: RECONSTRUCTION, COMPUTE MODEL!");
    return false;
  }

  const std::vector<Polygon_3>& planar_shapes() const {
    return m_polygons;
  }

  const Polygon_map& polygon_map() const {
    return m_polygon_map;
  }

  void clear() {
    m_ground_points.clear();
    m_boundary_points.clear();
    m_interior_points.clear();
    m_polygons.clear();
  }

private:
  const Input_range& m_input_range;
  const Semantic_map& m_semantic_map;

  Point_map_3  m_point_map_3;
  Vector_map_3 m_normal_map_3;

  Data_structure& m_data;
  const bool m_debug;
  const bool m_verbose;
  const Planar_shape_type m_planar_shape_type;
  const Converter m_converter;

  std::vector<std::size_t> m_ground_points;
  std::vector<std::size_t> m_boundary_points;
  std::vector<std::size_t> m_interior_points;

  std::vector<Polygon_3> m_polygons;
  Polygon_map m_polygon_map;

  std::map<std::size_t, Indices> m_region_map;

  void collect_points(
    const Semantic_label output_label,
    std::vector<std::size_t>& indices) const {

    indices.clear();
    for (std::size_t i = 0; i < m_input_range.size(); ++i) {
      const Semantic_label label =
      get(m_semantic_map, *(m_input_range.begin() + i));
      if (label == output_label)
        indices.push_back(i);
    }
  }

  void create_ground_plane() {

    if (m_verbose) std::cout << "* creating ground plane ... ";
    const auto plane = fit_plane(m_ground_points);
    const std::size_t shape_idx = add_planar_shape(m_ground_points, plane);
    CGAL_assertion(shape_idx != std::size_t(-1));
    m_region_map[shape_idx] = m_ground_points;
    if (m_verbose) std::cout << "done" << std::endl;
  }

  const Plane_3 fit_plane(const std::vector<std::size_t>& region) const {

    std::vector<IPoint_3> points;
    points.reserve(region.size());
    for (const std::size_t idx : region) {
      CGAL_assertion(idx < m_input_range.size());
      points.push_back(m_converter(get(m_point_map_3, idx)));
    }
    CGAL_assertion(points.size() == region.size());

    IPlane_3 fitted_plane;
    IPoint_3 fitted_centroid;
    CGAL::linear_least_squares_fitting_3(
      points.begin(), points.end(),
      fitted_plane, fitted_centroid,
      CGAL::Dimension_tag<0>());

    const Plane_3 plane(
      static_cast<FT>(fitted_plane.a()),
      static_cast<FT>(fitted_plane.b()),
      static_cast<FT>(fitted_plane.c()),
      static_cast<FT>(fitted_plane.d()));
    return plane;
  }

  const std::size_t add_planar_shape(
    const std::vector<std::size_t>& region, const Plane_3& plane) {

    switch (m_planar_shape_type) {
      case Planar_shape_type::CONVEX_HULL: {
        return add_convex_hull_shape(region, plane);
      }
      case Planar_shape_type::RECTANGLE: {
        return add_rectangle_shape(region, plane);
      }
      default: {
        CGAL_assertion_msg(false, "ERROR: ADD PLANAR SHAPE, WRONG TYPE!");
        return std::size_t(-1);
      }
    }
    return std::size_t(-1);
  }

  const std::size_t add_convex_hull_shape(
    const std::vector<std::size_t>& region, const Plane_3& plane) {

    std::vector<Point_2> points;
    points.reserve(region.size());
    for (const std::size_t idx : region) {
      CGAL_assertion(idx < m_input_range.size());
      const auto& p = get(m_point_map_3, idx);
      const auto q = plane.projection(p);
      const auto point = plane.to_2d(q);
      points.push_back(point);
    }
    CGAL_assertion(points.size() == region.size());

    std::vector<Point_2> ch;
    CGAL::convex_hull_2(points.begin(), points.end(), std::back_inserter(ch) );

    std::vector<Point_3> polygon;
    for (const auto& p : ch) {
      const auto point = plane.to_3d(p);
      polygon.push_back(point);
    }

    const std::size_t shape_idx = m_polygons.size();
    m_polygons.push_back(polygon);
    return shape_idx;
  }

  const std::size_t add_rectangle_shape(
    const std::vector<std::size_t>& region, const Plane_3& plane) {

    CGAL_assertion_msg(false, "TODO: ADD RECTANGLE SHAPE!");
    return std::size_t(-1);
  }

  template<typename NamedParameters>
  void create_approximate_walls(const NamedParameters& np) {

    std::vector< std::vector<std::size_t> > regions;
    apply_region_growing(np, m_boundary_points, regions);
    for (const auto& region : regions) {
      const auto plane = fit_plane(region);
      const std::size_t shape_idx = add_planar_shape(region, plane);
      CGAL_assertion(shape_idx != std::size_t(-1));
      m_region_map[shape_idx] = region;
    }
    if (m_verbose) {
      std::cout << "* found " << regions.size() << " approximate walls" << std::endl;
    }
  }

  template<typename NamedParameters>
  void create_approximate_roofs(const NamedParameters& np) {

    std::vector< std::vector<std::size_t> > regions;
    apply_region_growing(np, m_interior_points, regions);
    for (const auto& region : regions) {
      const auto plane = fit_plane(region);
      const std::size_t shape_idx = add_planar_shape(region, plane);
      CGAL_assertion(shape_idx != std::size_t(-1));
      m_region_map[shape_idx] = region;
    }
    if (m_verbose) {
      std::cout << "* found " << regions.size() << " approximate roofs" << std::endl;
    }
  }

  template<typename NamedParameters>
  void apply_region_growing(
    const NamedParameters& np,
    const std::vector<std::size_t>& input_range,
    std::vector< std::vector<std::size_t> >& regions) const {

    // const FT radius = parameters::choose_parameter(
    //   parameters::get_parameter(np, internal_np::neighbor_radius), FT(1));

    // Parameters.
    const std::size_t k = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::k_neighbors), 12);
    const FT max_distance_to_plane = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::distance_threshold), FT(1));
    const FT max_accepted_angle = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::angle_threshold), FT(15));
    const std::size_t min_region_size = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::min_region_size), 50);

    // Region growing.
    Neighbor_query_3 neighbor_query(input_range, k, m_point_map_3);

    Planar_region planar_region(input_range,
      max_distance_to_plane, max_accepted_angle, min_region_size,
      m_point_map_3, m_normal_map_3);

    Planar_sorting sorting(
      input_range, neighbor_query, m_point_map_3);
    sorting.sort();

    std::vector<Indices> result;
    Region_growing region_growing(
      input_range, neighbor_query, planar_region, sorting.seed_map());
    region_growing.detect(std::back_inserter(result));

    // Convert indices.
    regions.clear();
    regions.reserve(result.size());

    Indices region;
    for (const auto& indices : result) {
      region.clear();
      for (const std::size_t index : indices) {
        region.push_back(input_range[index]);
      }
      regions.push_back(region);
    }
    CGAL_assertion(regions.size() == result.size());
  }

  void assign_points_to_pfaces(std::map<PFace, Indices>& pface_points) const {

    CGAL_assertion(m_region_map.size() > 0);
    pface_points.clear();

    for (KSR::size_t i = 0; i < 6; ++i) {
      const auto pfaces = m_data.pfaces(i);
      for (const auto pface : pfaces) {
        pface_points[pface] = Indices();
      }
    }

    for (const auto& item : m_region_map) {
      const std::size_t shape_idx = item.first;
      const auto& indices = item.second;

      const KSR::size_t support_plane_idx = static_cast<KSR::size_t>(
        m_data.support_plane_index(shape_idx));
      CGAL_assertion(support_plane_idx >= 6);
      // dump_points(indices, "sp-points-" + std::to_string(support_plane_idx));

      const auto pfaces = m_data.pfaces(support_plane_idx);
      for (const auto pface : pfaces) {
        pface_points[pface] = Indices();
        const auto pvertices = m_data.pvertices_of_pface(pface);

        Delaunay tri;
        for (const auto pvertex : pvertices) {
          CGAL_assertion(m_data.has_ivertex(pvertex));
          const auto ivertex = m_data.ivertex(pvertex);
          const auto& point = m_data.point_2(support_plane_idx, ivertex);
          tri.insert(point);
        }

        for (const std::size_t index : indices) {
          const auto& point = get(m_point_map_3, index);
          const auto query = m_data.to_2d(support_plane_idx, point);
          const auto fh = tri.locate(query);
          if (fh != nullptr && !tri.is_infinite(fh)) {
            pface_points[pface].push_back(index);
          }
        }
      }
    }

    // for (const auto& item : pface_points) {
    //   dump_points(item.second, "pf-points-" + m_data.str(item.first));
    // }
  }

  void extract_surface() {
    CGAL_assertion_msg(false, "TODO: EXTRACT SURFACE FROM THE LABELED VOLUMES!");
  }

  void dump_points(
    const std::vector<std::size_t>& indices,
    const std::string file_name) const {

    std::vector<Point_3> points;
    points.reserve(indices.size());
    for (const std::size_t index : indices) {
      const auto& point = get(m_point_map_3, index);
      points.push_back(point);
    }
    CGAL_assertion(points.size() == indices.size());

    KSR_3::Saver<Kernel> saver;
    saver.export_points_3(points, file_name);
  }

  void dump_polygons(const std::string file_name) {

    KSR_3::Saver<Kernel> saver;
    saver.export_polygon_soup_3(m_polygons, file_name);
  }

  void dump_volumes(const std::string file_name) {

    for (const auto& volume : m_data.volumes()) {
      if (volume.visibility == Visibility_label::INSIDE) {
        dump_volume(m_data, volume.pfaces,
        file_name + "-" + std::to_string(volume.index), false);
      }
    }
  }
};

} // namespace KSR_3
} // namespace CGAL

#endif // CGAL_KSR_3_RECONSTRUCTION_H
