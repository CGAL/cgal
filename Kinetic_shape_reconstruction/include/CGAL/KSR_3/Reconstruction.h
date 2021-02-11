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
#include <CGAL/Regularization/regularize_planes.h>

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
  using FT       = typename Kernel::FT;
  using Point_2  = typename Kernel::Point_2;
  using Point_3  = typename Kernel::Point_3;
  using Plane_3  = typename Kernel::Plane_3;
  using Vector_3 = typename Kernel::Vector_3;

  using Data_structure = KSR_3::Data_structure<Kernel>;
  using PFace          = typename Data_structure::PFace;

  using Point_map_3        = KSR::Item_property_map<Input_range, Point_map>;
  using Vector_map_3       = KSR::Item_property_map<Input_range, Vector_map>;
  using Plane_map          = CGAL::Identity_property_map<Plane_3>;
  using Point_to_plane_map = CGAL::Shape_detection::RG::Point_to_shape_index_map;

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
  m_point_map(point_map),
  m_normal_map(normal_map),
  m_semantic_map(semantic_map),
  m_point_map_3(m_input_range, m_point_map),
  m_normal_map_3(m_input_range, m_normal_map),
  m_data(data),
  m_debug(true),
  m_verbose(verbose),
  m_planar_shape_type(Planar_shape_type::CONVEX_HULL) {

    clear();
    collect_points(Semantic_label::GROUND           , m_ground_points);
    collect_points(Semantic_label::BUILDING_BOUNDARY, m_boundary_points);
    collect_points(Semantic_label::BUILDING_INTERIOR, m_interior_points);

    // if (
    //   m_ground_points.size()   == 0 ||
    //   m_boundary_points.size() == 0 ||
    //   m_interior_points.size() == 0) {

    //   CGAL_assertion_msg(false, "TODO: IMPLEMENT FREE-FORM RECONSTRUCTION!");
    // }

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
    m_planes.clear();
    m_polygons.clear();
    m_region_map.clear();
    create_ground_plane();
    create_approximate_walls(np);
    create_approximate_roofs(np);
    CGAL_assertion(m_planes.size() == m_polygons.size());
    CGAL_assertion(m_polygons.size() == m_region_map.size());
    if (m_debug) dump_polygons("detected-planar-shapes");
    return true;
  }

  template<typename NamedParameters>
  const bool regularize_planar_shapes(
    const NamedParameters& np) {

    const bool regularize = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::regularize), false);
    if (!regularize) return true;

    // Regularize.
    const FT max_accepted_angle    = FT(10); // parameters::choose_parameter(
      // parameters::get_parameter(np, internal_np::angle_threshold), FT(15));
    const FT max_distance_to_plane = FT(1) / FT(5); // parameters::choose_parameter(
      // parameters::get_parameter(np, internal_np::distance_threshold), FT(1));
    const Vector_3 symmetry_axis(FT(0), FT(0), FT(1));

    if (m_verbose) {
      std::cout << std::endl << "--- REGULARIZING PLANAR SHAPES: " << std::endl;
    }

    std::vector<Plane_3> planes;
    std::vector<Indices> regions;
    create_planes_and_regions(planes, regions);

    CGAL_assertion(planes.size() > 0);
    CGAL_assertion(planes.size() == regions.size());

    Plane_map plane_map;
    Point_to_plane_map point_to_plane_map(m_input_range, regions);
    CGAL::regularize_planes(
      m_input_range,
      m_point_map,
      planes,
      plane_map,
      point_to_plane_map,
      true, true, true, false,
      max_accepted_angle,
      max_distance_to_plane,
      symmetry_axis);

    const std::size_t num_polygons = m_polygons.size();

    m_planes.clear();
    m_polygons.clear();
    m_region_map.clear();
    for (std::size_t i = 0; i < regions.size(); ++i) {
      const auto& plane  = planes[i];
      const auto& region = regions[i];

      const std::size_t shape_idx = add_planar_shape(region, plane);
      CGAL_assertion(shape_idx != std::size_t(-1));
      m_region_map[shape_idx] = region;
    }
    CGAL_assertion(m_polygons.size() == num_polygons);
    CGAL_assertion(m_polygons.size() == m_planes.size());
    CGAL_assertion(m_polygons.size() == m_region_map.size());

    if (m_verbose) {
      std::cout << "* num regularized planes: " << m_planes.size() << std::endl;
    }

    if (m_debug) dump_polygons("regularized-planar-shapes");
    return true;
  }

  template<typename NamedParameters>
  const bool compute_model(
    const NamedParameters& np) {

    if (m_verbose) {
      std::cout << std::endl << "--- COMPUTING THE MODEL: " << std::endl;
      std::cout << "* computing visibility ... ";
    }

    std::map<PFace, Indices> pface_points;
    assign_points_to_pfaces(pface_points);
    const Visibility visibility(
      m_data, pface_points, m_point_map_3, m_normal_map_3);

    CGAL_assertion(m_data.volumes().size() > 0);
    visibility.compute(m_data.volumes());
    // if (m_debug) dump_volumes("visibility/visibility");

    if (m_verbose) {
      std::cout << "done" << std::endl;
      std::cout << "* applying graphcut ... ";
    }

    const FT beta = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::graphcut_beta), FT(1) / FT(2));

    Graphcut graphcut(m_data, beta);
    graphcut.compute(m_data.volumes());
    // if (m_debug) dump_volumes("graphcut/graphcut");

    if (m_verbose) {
      std::cout << "done" << std::endl;
      std::cout << "* extracting the model ... ";
    }

    extract_surface_model();
    if (m_debug) dump_model("reconstructed-model");

    if (m_verbose) std::cout << "done" << std::endl;
    return true;
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
    m_planes.clear();
  }

private:
  const Input_range& m_input_range;
  const Point_map& m_point_map;
  const Vector_map& m_normal_map;
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
  std::vector<Plane_3> m_planes;
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

    if (m_ground_points.size() == 0) return;
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
    m_planes.push_back(plane);
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

  void create_planes_and_regions(
    std::vector<Plane_3>& planes,
    std::vector<Indices>& regions) const {

    planes.clear();
    planes.reserve(m_region_map.size());

    regions.clear();
    regions.reserve(m_region_map.size());

    CGAL_assertion(m_planes.size() == m_region_map.size());
    for (const auto& item : m_region_map) {
      const std::size_t shape_idx = item.first;

      const auto& plane = m_planes[shape_idx];
      CGAL_assertion(plane != Plane_3());
      planes.push_back(plane);

      const auto& region = item.second;
      CGAL_assertion(region.size() > 0);
      regions.push_back(region);
    }
    CGAL_assertion(planes.size()  == m_region_map.size());
    CGAL_assertion(regions.size() == m_region_map.size());
  }

  void assign_points_to_pfaces(std::map<PFace, Indices>& pface_points) const {

    pface_points.clear();
    for (std::size_t i = 0; i < m_data.number_of_support_planes(); ++i) {
      const auto pfaces = m_data.pfaces(i);
      for (const auto pface : pfaces) {
        pface_points[pface] = Indices();
      }
    }

    CGAL_assertion(m_region_map.size() > 0);
    for (const auto& item : m_region_map) {
      const std::size_t shape_idx = item.first;
      const auto& indices = item.second;

      const int sp_idx = m_data.support_plane_index(shape_idx);
      CGAL_assertion(sp_idx >= 6);
      const std::size_t support_plane_idx = static_cast<std::size_t>(sp_idx);
      // dump_points(indices, "sp-points-" + std::to_string(support_plane_idx));

      const auto pfaces = m_data.pfaces(support_plane_idx);
      for (const auto pface : pfaces) {
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

  void extract_surface_model() {
    create_surface_model();
    orient_surface_model();
  }

  void create_surface_model() {
    auto& model = m_data.reconstructed_model();
    model.clear();

    const auto& volumes = m_data.volumes();
    const auto& items   = m_data.pface_neighbors();

    for (const auto& item : items) {
      const auto& pface = item.first;
      const auto& neighbors = item.second;

      const int idx1 = neighbors.first;
      const int idx2 = neighbors.second;

      // std::cout << "idx1/2: " << idx1 << "/" << idx2 << std::endl;

      CGAL_assertion(idx1 >= 0 || idx2 >= 0);
      if (idx1 >= 0 && idx2 >= 0) {
        const auto& volume1 = volumes[idx1];
        const auto& volume2 = volumes[idx2];

        const auto label1 = volume1.visibility;
        const auto label2 = volume2.visibility;

        if (
          label1 == Visibility_label::INSIDE &&
          label2 == Visibility_label::OUTSIDE) {
          model.pfaces.push_back(pface);
        } else if (
          label1 == Visibility_label::OUTSIDE &&
          label2 == Visibility_label::INSIDE) {
          model.pfaces.push_back(pface);
        }
        continue;
      }

      if (idx1 >= 0) {
        CGAL_assertion(idx2 < 0);
        const auto& volume1 = volumes[idx1];
        const auto label1 = volume1.visibility;
        if (label1 == Visibility_label::INSIDE) {
          model.pfaces.push_back(pface);
        }
        continue;
      }

      if (idx2 >= 0) {
        CGAL_assertion(idx1 < 0);
        const auto& volume2 = volumes[idx2];
        const auto label2 = volume2.visibility;
        if (label2 == Visibility_label::INSIDE) {
          model.pfaces.push_back(pface);
        }
        continue;
      }
    }
  }

  void orient_surface_model() {
    return;
    CGAL_assertion_msg(false, "TODO: ORIENT SURFACE MODEL!");
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

  void dump_model(const std::string file_name) {

    std::vector<Point_3> polygon;
    std::vector< std::vector<Point_3> > polygons;
    const auto& model = m_data.reconstructed_model();
    std::vector<CGAL::Color> colors;

    std::size_t polygon_id = 0;
    KSR_3::Saver<Kernel> saver;
    for (const auto& pface : model.pfaces) {
      const auto pvertices = m_data.pvertices_of_pface(pface);
      polygon.clear();
      for (const auto pvertex : pvertices) {
        CGAL_assertion(m_data.has_ivertex(pvertex));
        const auto ivertex = m_data.ivertex(pvertex);
        const auto point = m_data.point_3(ivertex);
        polygon.push_back(point);
      }
      polygons.push_back(polygon);
      colors.push_back(saver.get_idx_color(pface.first));
      ++polygon_id;
    }
    saver.export_polygon_soup_3(polygons, colors, file_name);
  }
};

} // namespace KSR_3
} // namespace CGAL

#endif // CGAL_KSR_3_RECONSTRUCTION_H
