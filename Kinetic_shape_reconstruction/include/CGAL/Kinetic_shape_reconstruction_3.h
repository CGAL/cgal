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

#ifndef CGAL_KINETIC_SHAPE_RECONSTRUCTION_3_H
#define CGAL_KINETIC_SHAPE_RECONSTRUCTION_3_H

// #include <CGAL/license/Kinetic_shape_reconstruction.h>

// Boost includes.
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/boost/graph/Named_function_parameters.h>

// CGAL includes.
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>

// Internal includes.
#include <CGAL/KSR/utils.h>
#include <CGAL/KSR/debug.h>
#include <CGAL/KSR_3/Event.h>
#include <CGAL/KSR_3/Event_queue.h>
#include <CGAL/KSR_3/Data_structure.h>
#include <CGAL/KSR_3/Reconstruction.h>
#include <CGAL/KSR_3/Initializer.h>

namespace CGAL {

template<typename GeomTraits>
class Kinetic_shape_reconstruction_3 {

public:
  using Kernel = GeomTraits;

private:
  using FT        = typename Kernel::FT;
  using Point_2   = typename Kernel::Point_2;
  using Point_3   = typename Kernel::Point_3;
  using Vector_2  = typename Kernel::Vector_2;
  using Segment_2 = typename Kernel::Segment_2;

  using Data_structure = KSR_3::Data_structure<Kernel>;

  using PVertex = typename Data_structure::PVertex;
  using PFace   = typename Data_structure::PFace;

  using IVertex = typename Data_structure::IVertex;
  using IEdge   = typename Data_structure::IEdge;

  using Event       = KSR_3::Event<Data_structure>;
  using Event_queue = KSR_3::Event_queue<Data_structure>;

  using Bbox_2      = CGAL::Bbox_2;
  using EK          = CGAL::Exact_predicates_exact_constructions_kernel;
  using Initializer = KSR_3::Initializer<EK>;

  using Polygon_mesh = CGAL::Surface_mesh<Point_3>;
  using Vertex_index = typename Polygon_mesh::Vertex_index;

private:
  const bool m_debug;
  const bool m_verbose;
  Event_queue m_queue;
  FT m_min_time;
  FT m_max_time;
  Initializer m_initializer;
  Data_structure m_data;

public:
  Kinetic_shape_reconstruction_3(
    const bool verbose = true,
    const bool debug   = false) :
  m_debug(debug),
  m_verbose(verbose),
  m_queue(m_debug),
  m_min_time(-FT(1)),
  m_max_time(-FT(1)),
  m_initializer(m_debug, m_verbose),
  m_data(m_debug)
  { }

  template<
  typename InputRange,
  typename PolygonMap,
  typename NamedParameters>
  const bool partition(
    const InputRange& input_range,
    const PolygonMap polygon_map,
    const NamedParameters& np) {

    const unsigned int k = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::k_intersections), 1);
    unsigned int n = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::n_subdivisions), 0);
    double enlarge_bbox_ratio = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::enlarge_bbox_ratio), 1.1);
    const bool reorient = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::reorient), false);

    std::cout.precision(20);
    if (input_range.size() == 0) {
      CGAL_warning_msg(input_range.size() > 0,
      "WARNING: YOUR INPUT IS EMPTY! RETURN WITH NO CHANGE!");
      return false;
    }

    if (n != 0) {
      CGAL_assertion_msg(false, "TODO: IMPLEMENT KINETIC SUBDIVISION!");
      if (n > 3) {
        CGAL_warning_msg(n <= 3,
        "WARNING: DOES IT MAKE SENSE TO HAVE MORE THAN 64 INPUT BLOCKS? SETTING N TO 3!");
        n = 3;
      }
    }

    if (enlarge_bbox_ratio < 1.0) {
      CGAL_warning_msg(enlarge_bbox_ratio >= 1.0,
      "WARNING: YOU SET ENLARGE_BBOX_RATIO < 1.0! THE VALID RANGE IS [1.0, +INF). SETTING TO 1.0!");
      enlarge_bbox_ratio = 1.0;
    }

    if (m_verbose) {
      const unsigned int num_blocks = std::pow(n + 1, 3);
      const std::string is_reorient = (reorient ? "true" : "false");

      std::cout << std::endl << "--- PARTITION OPTIONS: " << std::endl;
      std::cout << "* number of intersections k: "            << k                  << std::endl;
      std::cout << "* number of subdivisions per bbox side: " << n                  << std::endl;
      std::cout << "* number of subdivision blocks: "         << num_blocks         << std::endl;
      std::cout << "* enlarge bbox ratio: "                   << enlarge_bbox_ratio << std::endl;
      std::cout << "* reorient: "                             << is_reorient        << std::endl;
    }

    const FT time_step = static_cast<FT>(m_initializer.initialize(
      input_range, polygon_map, k, enlarge_bbox_ratio, reorient));
    m_initializer.convert(m_data);
    m_data.check_integrity();

    if (k == 0) {
      CGAL_warning_msg(k > 0,
      "WARNING: YOU SET K TO 0! THAT MEANS NO PROPAGATION! THE VALID VALUES ARE {1,2,...}. INTERSECT AND RETURN!");
      return false;
    }

    // if (m_verbose) {
    //   std::cout << std::endl << "POLYGON SPLITTER SUCCESS!" << std::endl << std::endl;
    // }
    // return true;
    // exit(EXIT_SUCCESS);

    if (m_verbose) {
      std::cout << std::endl << "--- RUNNING THE QUEUE:" << std::endl;
      std::cout << "* propagation started" << std::endl;
    }
    std::size_t num_iterations = 0;
    m_min_time = FT(0);
    m_max_time = time_step;
    CGAL_assertion(m_min_time >= FT(0) && m_max_time >= m_min_time);
    std::size_t global_iteration = 0;
    while (initialize_queue()) {

      global_iteration = run(k, global_iteration);
      m_min_time = m_max_time;
      m_max_time += time_step;
      m_data.check_integrity();
      ++num_iterations;

      if (m_verbose && !m_debug) {
        std::cout << ".";
        if (num_iterations == 50) {
          std::cout << std::endl;
        }
      }

      if (num_iterations > 100000000) {
        CGAL_assertion_msg(false, "DEBUG WARNING: WHY SO MANY ITERATIONS?");
      }
    }
    if (m_verbose) {
      if (m_verbose && !m_debug) std::cout << std::endl;
      std::cout << "* propagation finished" << std::endl;
      std::cout << "* number of events: " << global_iteration << std::endl;
    }

    if (m_verbose) std::cout << std::endl << "--- FINALIZING PARTITION:" << std::endl;
    if (m_debug) dump(m_data, "jiter-final-a-result");
    m_data.finalize();
    if (m_verbose) std::cout << "* checking final mesh integrity ...";
    m_data.check_integrity();
    if (m_verbose) std::cout << " done" << std::endl;
    if (m_debug) dump(m_data, "jiter-final-b-result");

    // std::cout << std::endl << "CLEANING SUCCESS!" << std::endl << std::endl;
    // exit(EXIT_SUCCESS);

    if (m_verbose) std::cout << "* getting volumes:" << std::endl;
    m_data.create_polyhedra();
    return true;
  }

  template<
  typename InputRange,
  typename PointMap,
  typename VectorMap,
  typename SemanticMap,
  typename NamedParameters>
  const bool reconstruct(
    const InputRange& input_range,
    const PointMap point_map,
    const VectorMap normal_map,
    const SemanticMap semantic_map,
    const NamedParameters& np) {

    using Reconstruction = KSR_3::Reconstruction<
      InputRange, PointMap, VectorMap, SemanticMap, Kernel>;

    Reconstruction reconstruction(
      input_range, point_map, normal_map, semantic_map, m_data, m_verbose);
    bool success = reconstruction.detect_planar_shapes(np);
    if (!success) {
      CGAL_assertion_msg(false, "ERROR: RECONSTRUCTION, DETECTING PLANAR SHAPES FAILED!");
    }

    success = partition(
      reconstruction.planar_shapes(), reconstruction.polygon_map(), np);
    if (!success) {
      CGAL_assertion_msg(false, "ERROR: RECONSTRUCTION, PARTITION FAILED!");
    }

    success = reconstruction.compute_model(np);
    if (!success) {
      CGAL_assertion_msg(false, "ERROR: RECONSTRUCTION, COMPUTING MODEL FAILED!");
    }
    return success;
  }

  /*******************************
  **         STATISTICS         **
  ********************************/

  const int number_of_support_planes() const {
    return static_cast<int>(m_data.number_of_support_planes());
  }

  const std::size_t number_of_vertices(const int support_plane_idx = -1) const {

    CGAL_assertion(support_plane_idx < number_of_support_planes());
    if (support_plane_idx >= number_of_support_planes()) return std::size_t(-1);
    if (support_plane_idx < 0) {
      return m_data.igraph().number_of_vertices();
    }

    CGAL_assertion(support_plane_idx >= 0);
    const KSR::size_t sp_idx = static_cast<KSR::size_t>(support_plane_idx);
    return static_cast<std::size_t>(m_data.mesh(sp_idx).number_of_vertices());
  }

  const std::size_t number_of_edges(const int support_plane_idx = -1) const {

    CGAL_assertion(support_plane_idx < number_of_support_planes());
    if (support_plane_idx >= number_of_support_planes()) return std::size_t(-1);
    if (support_plane_idx < 0) {
      return m_data.igraph().number_of_edges();
    }

    CGAL_assertion(support_plane_idx >= 0);
    const KSR::size_t sp_idx = static_cast<KSR::size_t>(support_plane_idx);
    return static_cast<std::size_t>(m_data.mesh(sp_idx).number_of_edges());
  }

  const std::size_t number_of_faces(const int support_plane_idx = -1) const {

    CGAL_assertion(support_plane_idx < number_of_support_planes());
    if (support_plane_idx >= number_of_support_planes()) return std::size_t(-1);
    if (support_plane_idx < 0) {
      std::size_t num_all_faces = 0;
      for (int i = 0; i < number_of_support_planes(); ++i) {
        const std::size_t num_faces = static_cast<std::size_t>(
          m_data.mesh(static_cast<KSR::size_t>(i)).number_of_faces());
        num_all_faces += num_faces;
      }
      return num_all_faces;
    }

    CGAL_assertion(support_plane_idx >= 0);
    const KSR::size_t sp_idx = static_cast<KSR::size_t>(support_plane_idx);
    const std::size_t num_faces = static_cast<std::size_t>(
      m_data.mesh(sp_idx).number_of_faces());
    return num_faces;
  }

  const int number_of_volume_levels() const {
    return m_data.number_of_volume_levels();
  }

  const std::size_t number_of_volumes(const int volume_level = -1) const {
    return m_data.number_of_volumes(volume_level);
  }

  const int support_plane_index(const std::size_t polygon_index) const {
    const int support_plane_idx = m_data.support_plane_index(polygon_index);
    CGAL_assertion(support_plane_idx >= 6);
    return support_plane_idx;
  }

  /*******************************
  **           OUTPUT           **
  ********************************/

  template<typename VertexOutputIterator>
  VertexOutputIterator output_partition_vertices(
    VertexOutputIterator vertices, const int support_plane_idx = -1) const {

    CGAL_assertion(support_plane_idx < number_of_support_planes());
    if (support_plane_idx >= number_of_support_planes()) return vertices;
    if (support_plane_idx < 0) {
      const auto all_ivertices = m_data.ivertices();
      for (const auto ivertex : all_ivertices) {
        *(vertices++) = m_data.point_3(ivertex);
      }
      return vertices;
    }

    CGAL_assertion(support_plane_idx >= 0);
    const KSR::size_t sp_idx = static_cast<KSR::size_t>(support_plane_idx);
    const auto all_pvertices = m_data.pvertices(sp_idx);
    for (const auto pvertex : all_pvertices) {
      CGAL_assertion(m_data.has_ivertex(pvertex));
      const auto ivertex = m_data.ivertex(pvertex);
      *(vertices++) = m_data.point_3(ivertex);
    }
    return vertices;
  }

  template<typename EdgeOutputIterator>
  EdgeOutputIterator output_partition_edges(
    EdgeOutputIterator edges, const int support_plane_idx = -1) const {

    CGAL_assertion(support_plane_idx < number_of_support_planes());
    if (support_plane_idx >= number_of_support_planes()) return edges;
    if (support_plane_idx < 0) {
      const auto all_iedges = m_data.iedges();
      for (const auto iedge : all_iedges) {
        *(edges++) = m_data.segment_3(iedge);
      }
      return edges;
    }

    CGAL_assertion(support_plane_idx >= 0);
    const KSR::size_t sp_idx = static_cast<KSR::size_t>(support_plane_idx);
    const auto all_pedges = m_data.pedges(sp_idx);
    for (const auto pedge : all_pedges) {
      CGAL_assertion(m_data.has_iedge(pedge));
      const auto iedge = m_data.iedge(pedge);
      *(edges++) = m_data.segment_3(iedge);
    }
    return edges;
  }

  template<typename FaceOutputIterator>
  FaceOutputIterator output_partition_faces(
    FaceOutputIterator faces, const int support_plane_idx = -1) const {

    KSR::Indexer<IVertex> indexer;
    CGAL_assertion(support_plane_idx < number_of_support_planes());
    if (support_plane_idx >= number_of_support_planes()) return faces;
    if (support_plane_idx < 0) {
      const auto all_ivertices = m_data.ivertices();
      for (const auto ivertex : all_ivertices) indexer(ivertex);
      for (int i = 0; i < number_of_support_planes(); ++i) {
        const KSR::size_t sp_idx = static_cast<KSR::size_t>(i);
        output_partition_faces(faces, indexer, sp_idx);
      }
      return faces;
    }

    CGAL_assertion(support_plane_idx >= 0);
    const KSR::size_t sp_idx = static_cast<KSR::size_t>(support_plane_idx);
    const auto all_pvertices = m_data.pvertices(sp_idx);
    for (const auto pvertex : all_pvertices) {
      CGAL_assertion(m_data.has_ivertex(pvertex));
      const auto ivertex = m_data.ivertex(pvertex);
      indexer(ivertex);
    }
    return output_partition_faces(faces, indexer, sp_idx);
  }

  void output_support_plane(
    Polygon_mesh& polygon_mesh, const int support_plane_idx) const {

    polygon_mesh.clear();
    CGAL_assertion(support_plane_idx >= 0);
    if (support_plane_idx < 0) return;
    CGAL_assertion(support_plane_idx < number_of_support_planes());
    if (support_plane_idx >= number_of_support_planes()) return;
    const KSR::size_t sp_idx = static_cast<KSR::size_t>(support_plane_idx);

    std::vector<Vertex_index> vertices;
    std::vector<Vertex_index> map_vertices;

    map_vertices.clear();
    const auto all_pvertices = m_data.pvertices(sp_idx);
    for (const auto pvertex : all_pvertices) {
      CGAL_assertion(m_data.has_ivertex(pvertex));
      const auto ivertex = m_data.ivertex(pvertex);

      if (map_vertices.size() <= pvertex.second)
        map_vertices.resize(pvertex.second + 1);
      map_vertices[pvertex.second] =
        polygon_mesh.add_vertex(m_data.point_3(ivertex));
    }

    const auto all_pfaces = m_data.pfaces(sp_idx);
    for (const auto pface : all_pfaces) {
      vertices.clear();
      const auto pvertices = m_data.pvertices_of_pface(pface);
      for (const auto pvertex : pvertices) {
        vertices.push_back(map_vertices[pvertex.second]);
      }
      polygon_mesh.add_face(vertices);
    }
  }

  template<typename VolumeOutputIterator>
  VolumeOutputIterator output_partition_volumes(
    VolumeOutputIterator volumes, const int volume_level = -1) const {

    CGAL_assertion(volume_level < number_of_volume_levels());
    if (volume_level >= number_of_volume_levels()) return volumes;
    if (volume_level < 0) {
      for (int i = 0; i < number_of_volume_levels(); ++i) {
        output_partition_volumes(volumes, i);
      }
      return volumes;
    }

    CGAL_assertion(volume_level >= 0);
    std::size_t begin = 0;
    if (volume_level > 0) {
      for (int i = 0; i < volume_level; ++i) {
        begin += number_of_volumes(i);
      }
    }
    const std::size_t end = begin + number_of_volumes(volume_level);
    for (std::size_t i = begin; i < end; ++i) {
      output_partition_volume(volumes, i);
    }
    return volumes;
  }

  template<typename VolumeOutputIterator>
  VolumeOutputIterator output_partition_volume(
    VolumeOutputIterator volumes, const std::size_t volume_index) const {

    CGAL_assertion(volume_index < number_of_volumes(-1));
    if (volume_index >= number_of_volumes(-1)) return volumes;

    std::vector<Point_3> vertices;
    std::vector< std::vector<std::size_t> > faces;
    output_partition_volume(
      std::back_inserter(vertices), std::back_inserter(faces), volume_index);
    CGAL::Polygon_mesh_processing::orient_polygon_soup(vertices, faces);

    Polygon_mesh polygon_mesh;
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(
      vertices, faces, polygon_mesh);
    *(volumes++) = polygon_mesh;
    return volumes;
  }

  template<typename VertexOutputIterator, typename FaceOutputIterator>
  void output_partition_volume(
    VertexOutputIterator vertices, FaceOutputIterator faces,
    const std::size_t volume_index) const {

    CGAL_assertion(volume_index < number_of_volumes(-1));
    if (volume_index >= number_of_volumes(-1)) return;
    CGAL_assertion(m_data.volumes().size() == number_of_volumes(-1));
    const auto& volume = m_data.volumes()[volume_index];

    std::size_t num_vertices = 0;
    KSR::Indexer<IVertex> indexer;

    std::vector<std::size_t> face;
    const auto& pfaces = volume.pfaces;
    for (const auto& pface : pfaces) {
      face.clear();
      const auto pvertices = m_data.pvertices_of_pface(pface);
      for (const auto pvertex : pvertices) {

        CGAL_assertion(m_data.has_ivertex(pvertex));
        const auto ivertex = m_data.ivertex(pvertex);
        const std::size_t idx = indexer(ivertex);

        if (idx == num_vertices) {
          *(vertices++) = m_data.point_3(ivertex);
          ++num_vertices;
        }
        face.push_back(idx);
      }
      *(faces++) = face;
    }
  }

  template<typename LCC>
  void output_partition(LCC& lcc) const {

    CGAL_assertion_msg(false, "TODO: OUTPUT PARTITION LCC!");
  }

  template<typename VertexOutputIterator, typename FaceOutputIterator>
  void output_reconstructed_model(
    VertexOutputIterator vertices, FaceOutputIterator faces) const {

    CGAL_assertion_msg(false, "TODO: OUTPUT RECONSTRUCTED MODEL!");
  }

  /*******************************
  **           MEMORY           **
  ********************************/

  void clear() {
    m_data.clear();
    m_queue.clear();
    m_min_time = -FT(1);
    m_max_time = -FT(1);
    m_initializer.clear();
  }

private:

  template<typename FaceOutputIterator>
  FaceOutputIterator output_partition_faces(
    FaceOutputIterator faces, KSR::Indexer<IVertex>& indexer,
    const KSR::size_t sp_idx) const {

    std::vector<std::size_t> face;
    const auto all_pfaces = m_data.pfaces(sp_idx);
    for (const auto pface : all_pfaces) {
      face.clear();
      const auto pvertices = m_data.pvertices_of_pface(pface);
      for (const auto pvertex : pvertices) {
        CGAL_assertion(m_data.has_ivertex(pvertex));
        const auto ivertex = m_data.ivertex(pvertex);
        const std::size_t idx = indexer(ivertex);
        face.push_back(idx);
      }
      *(faces++) = face;
    }
    return faces;
  }

  const bool initialize_queue() {

    if (m_debug) {
      std::cout << "* initializing queue for events in [" <<
      m_min_time << ";" << m_max_time << "]" << std::endl;
    }

    m_data.update_positions(m_max_time);
    bool still_running = false;

    KSR::vector<IEdge> iedges;
    KSR::vector<Segment_2> segments;
    KSR::vector<Bbox_2> bboxes;
    for (KSR::size_t i = 0; i < m_data.number_of_support_planes(); ++i) {
      initialize_search_structures(i, iedges, segments, bboxes);
      for (const auto pvertex : m_data.pvertices(i)) {
        if (compute_events_of_pvertex(pvertex, iedges, segments, bboxes)) {
          still_running = true;
        }
      }
    }
    m_data.update_positions(m_min_time);
    return still_running;
  }

  void initialize_search_structures(
    const KSR::size_t i,
    KSR::vector<IEdge>& iedges,
    KSR::vector<Segment_2>& segments,
    KSR::vector<Bbox_2>& bboxes) {

    iedges.clear();
    segments.clear();
    bboxes.clear();

    // To get random access, copy in vector (suboptimal to do this
    // all the time, maybe this should be done once and for all and
    // replace the set).
    iedges.reserve(m_data.iedges(i).size());
    std::copy(m_data.iedges(i).begin(), m_data.iedges(i).end(), std::back_inserter(iedges));

    // Precompute segments and bboxes.
    segments.reserve(iedges.size());
    bboxes.reserve(iedges.size());
    for (const auto& iedge : iedges) {
      segments.push_back(m_data.segment_2(i, iedge));
      bboxes.push_back(segments.back().bbox());
    }
  }

  template<typename PVertexRange>
  void compute_events_of_pvertices(
    const FT last_event_time, const PVertexRange& pvertices) {

    m_min_time = m_data.current_time();
    m_data.update_positions(m_max_time);

    KSR::vector<IEdge> iedges;
    KSR::vector<Segment_2> segments;
    KSR::vector<Bbox_2> bboxes;
    initialize_search_structures(
      pvertices.front().first, iedges, segments, bboxes);

    for (const auto& pvertex : pvertices) {
      m_data.deactivate(pvertex);
    }

    for (const auto& pvertex : pvertices) {
      m_data.set_last_event_time(pvertex, last_event_time);
      compute_events_of_pvertex(pvertex, iedges, segments, bboxes);
    }

    for (const auto& pvertex : pvertices) {
      m_data.activate(pvertex);
    }

    m_data.update_positions(m_min_time);
  }

  const bool compute_events_of_pvertex(
    const PVertex& pvertex,
    const KSR::vector<IEdge>& iedges,
    const KSR::vector<Segment_2>& segments,
    const KSR::vector<Bbox_2>& bboxes) {

    std::cout.precision(20);
    if (m_data.is_frozen(pvertex)) {
      return false;
    }

    const auto pv_min = m_data.point_2(pvertex, m_min_time);
    const auto pv_max = m_data.point_2(pvertex, m_max_time);
    const Segment_2 pv_segment(pv_min, pv_max);
    const auto pv_bbox = pv_segment.bbox();

    if (m_data.has_iedge(pvertex)) {
      compute_events_of_constrained_pvertex(
        pvertex, pv_segment, pv_bbox);
    } else {
      compute_events_of_unconstrained_pvertex(
        pvertex, pv_segment, pv_bbox, iedges, segments, bboxes);
    }
    return true;
  }

  void compute_events_of_constrained_pvertex(
    const PVertex& pvertex,
    const Segment_2& pv_segment,
    const Bbox_2& pv_bbox) {

    Event event1, event2;
    bool is_pvertex_to_pvertex_event, is_pvertex_to_ivertex_event;
    std::tie(is_pvertex_to_pvertex_event, event1) =
      try_pvertex_to_pvertex_constrained_event(pvertex, pv_segment, pv_bbox);
    std::tie(is_pvertex_to_ivertex_event, event2) =
      try_pvertex_to_ivertex_constrained_event(pvertex, pv_segment, pv_bbox);

    if (is_pvertex_to_pvertex_event) {
      m_queue.push(event1);
    }
    if (is_pvertex_to_ivertex_event) {
      m_queue.push(event2);
    }

    // Is this version better? Does it work?
    // if (is_pvertex_to_pvertex_event && is_pvertex_to_ivertex_event) {
    //   if (event1.time() < event2.time()) {
    //     m_queue.push(event1);
    //   } else if (event2.time() < event1.time()) {
    //     m_queue.push(event2);
    //   } else {
    //     CGAL_assertion_msg(false, "TODO: WHAT SHOULD WE DO FOR EQUAL EVENTS?");
    //   } return;
    // }
    // if (is_pvertex_to_pvertex_event) {
    //   CGAL_assertion(!is_pvertex_to_ivertex_event);
    //   m_queue.push(event1); return;
    // }
    // if (is_pvertex_to_ivertex_event) {
    //   CGAL_assertion(!is_pvertex_to_pvertex_event);
    //   m_queue.push(event2); return;
    // }
    // CGAL_assertion(!is_pvertex_to_pvertex_event && !is_pvertex_to_ivertex_event);
  }

  // Test left and right vertices of the mesh face.
  const std::pair<bool, Event> try_pvertex_to_pvertex_constrained_event(
    const PVertex& pvertex,
    const Segment_2& pv_segment,
    const Bbox_2& pv_bbox) {

    Event event;
    bool is_event_found = false;

    PVertex prev, next;
    std::tie(prev, next) = m_data.prev_and_next(pvertex);
    for (const auto& pother : { prev, next }) {
      if (pother == Data_structure::null_pvertex()
          || !m_data.is_active(pother)
          ||  m_data.has_iedge(pother)) {
        continue;
      }

      const Segment_2 po_segment(
        m_data.point_2(pother, m_min_time),
        m_data.point_2(pother, m_max_time));
      const auto po_bbox = po_segment.bbox();

      if (!do_overlap(pv_bbox, po_bbox)) {
        continue;
      }

      Point_2 inter;
      if (!KSR::intersection(pv_segment, po_segment, inter)) {
        continue;
      }

      const FT distance = KSR::distance(pv_segment.source(), inter);
      const FT time = distance / m_data.speed(pvertex);

      // Constrained pvertex to another pvertex event.
      event = Event(true, pvertex, pother, m_min_time + time);
      is_event_found = true;

      // std::cout << "pvertex: " << m_data.point_3(pvertex) << std::endl;
      // std::cout << "pother: "  << m_data.point_3(pother)  << std::endl;
    }
    return std::make_pair(is_event_found, event);
  }

  // Test end vertices of the intersection edge.
  const std::pair<bool, Event> try_pvertex_to_ivertex_constrained_event(
    const PVertex& pvertex,
    const Segment_2& pv_segment,
    const Bbox_2& pv_bbox) {

    Event event;
    bool is_event_found = false;

    CGAL_assertion(m_data.has_iedge(pvertex));
    const auto iedge = m_data.iedge(pvertex);
    for (const auto& ivertex : { m_data.source(iedge), m_data.target(iedge) }) {
      if (!m_data.is_active(ivertex)) {
        continue;
      }

      const Point_2 ipoint = m_data.to_2d(pvertex.first, ivertex);
      const auto vec1 = pv_segment.to_vector();
      const Vector_2 vec2(pv_segment.source(), ipoint);
      const FT dot_product = vec1 * vec2;
      if (dot_product < FT(0)) { // opposite directions
        continue;
      }

      const FT distance = KSR::distance(pv_segment.source(), ipoint);
      const FT time = distance / m_data.speed(pvertex);

      // Constrained pvertex to ivertex event.
      if (time < m_max_time - m_min_time) {
        event = Event(true, pvertex, ivertex, m_min_time + time);
        is_event_found = true;

        // std::cout << "pvertex: " << m_data.point_3(pvertex) << std::endl;
        // std::cout << "ivertex: " << m_data.point_3(ivertex) << std::endl;
      }
    }
    return std::make_pair(is_event_found, event);
  }

  void compute_events_of_unconstrained_pvertex(
    const PVertex& pvertex,
    const Segment_2& pv_segment,
    const Bbox_2& pv_bbox,
    const KSR::vector<IEdge>& iedges,
    const KSR::vector<Segment_2>& segments,
    const KSR::vector<Bbox_2>& bboxes) {

    try_pvertex_to_iedge_unconstrained_event(
      pvertex, pv_segment, pv_bbox, iedges, segments, bboxes);
  }

  // Test all intersection edges.
  const bool try_pvertex_to_iedge_unconstrained_event(
    const PVertex& pvertex,
    const Segment_2& pv_segment,
    const Bbox_2& pv_bbox,
    const KSR::vector<IEdge>& iedges,
    const KSR::vector<Segment_2>& segments,
    const KSR::vector<Bbox_2>& bboxes) {

    bool is_event_found = false;
    const auto prev = m_data.prev(pvertex);
    const auto next = m_data.next(pvertex);
    for (std::size_t i = 0; i < iedges.size(); ++i) {
      const auto& iedge = iedges[i];

      if (m_data.iedge(prev) == iedge ||
          m_data.iedge(next) == iedge) {
        continue;
      }

      if (!m_data.is_active(iedge)) {
        continue;
      }

      if (!CGAL::do_overlap(pv_bbox, bboxes[i])) {
        continue;
      }

      Point_2 inter;
      if (!KSR::intersection(pv_segment, segments[i], inter)) {
        continue;
      }

      // Try to add unconstrained pvertex to ivertex event.
      is_event_found = try_pvertex_to_ivertex_unconstrained_event(
        pvertex, iedge, inter);

      // Otherwise we add unconstrained pvertex to iedge event.
      if (!is_event_found) {
        const FT distance = KSR::distance(pv_segment.source(), inter);
        const FT time = distance / m_data.speed(pvertex);
        m_queue.push(Event(false, pvertex, iedge, m_min_time + time));
        is_event_found = true;
      }

      // std::cout << "pvertex: " << m_data.point_3(pvertex) << std::endl;
      // std::cout << "iedge: "   << m_data.segment_3(iedge) << std::endl;
    }
    return is_event_found;
  }

  const bool try_pvertex_to_ivertex_unconstrained_event(
    const PVertex& pvertex,
    const IEdge& iedge,
    const Point_2& inter) {

    const bool is_event_found = false;
    return is_event_found;
    CGAL_assertion_msg(false,
    "TODO: ADD PVERTEX TO IVERTEX UNCONSTRAINED EVENT!");
  }

  const std::size_t run(
    const unsigned int k,
    const std::size_t initial_iteration) {

    if (m_debug) {
      std::cout << "* unstacking queue, current size: " << m_queue.size() << std::endl;
    }

    std::size_t iteration = initial_iteration;
    while (!m_queue.empty()) {

      const Event event = m_queue.pop();
      const FT current_time = event.time();
      if (m_debug) {
        if (iteration < 10) {
          dump(m_data, "iter-0" + std::to_string(iteration));
          dump_event(m_data, event, "iter-0" + std::to_string(iteration));
        } else {
          dump(m_data, "iter-" + std::to_string(iteration));
          dump_event(m_data, event, "iter-" + std::to_string(iteration));
        }
        // const KSR::size_t sp_debug_idx = 23;
        // dump_2d_surface_mesh(m_data, sp_debug_idx, "iter-" + std::to_string(iteration) +
        // "-surface-mesh-" + std::to_string(sp_debug_idx));
      }

      m_data.update_positions(current_time);
      if (m_debug) {
        std::cout << std::endl << "* APPLYING " << iteration << ": " << event << std::endl;
      }
      ++iteration;

      // if (iteration == 380) {
      //   exit(EXIT_FAILURE);
      // }

      apply(event, k);
      m_data.check_integrity();
    }
    return iteration;
  }

  void apply(const Event& event, const unsigned int k) {

    const auto pvertex = event.pvertex();
    if (event.is_pvertex_to_pvertex()) {
      const auto pother = event.pother();

      remove_events(pvertex);
      remove_events(pother);

      if (m_data.has_iedge(pvertex)) {
        CGAL_assertion(m_data.has_iedge(pvertex));
        if (m_data.has_iedge(pother)) {
          apply_event_two_constrained_pvertices_meet(pvertex, pother, event);
        } else {
          apply_event_constrained_pvertex_meets_free_pvertex(pvertex, pother, event);
        }
      } else {
        CGAL_assertion(!m_data.has_iedge(pvertex));
        if (!m_data.has_iedge(pother)) {
          apply_event_two_unconstrained_pvertices_meet(pvertex, pother, event);
        } else {
          CGAL_assertion_msg(false, "ERROR: THIS EVENT SHOULD NOT EVER HAPPEN!");
          apply_event_constrained_pvertex_meets_free_pvertex(pother, pvertex, event);
        }
      }
    } else if (event.is_pvertex_to_iedge()) {

      const auto iedge = event.iedge();
      if (m_data.has_iedge(pvertex)) {
        apply_event_constrained_pvertex_meets_iedge(pvertex, iedge, event);
      } else {
        const bool is_event_happend = apply_event_unconstrained_pedge_meets_iedge(
          pvertex, iedge, event);
        if (!is_event_happend) {
          apply_event_unconstrained_pvertex_meets_iedge(pvertex, iedge, event);
        }
      }
    } else if (event.is_pvertex_to_ivertex()) {

      const auto ivertex = event.ivertex();
      if (m_data.has_iedge(pvertex)) {
        apply_event_constrained_pvertex_meets_ivertex(pvertex, ivertex, event);
      } else {
        apply_event_unconstrained_pvertex_meets_ivertex(pvertex, ivertex, event);
      }
    } else {
      CGAL_assertion_msg(false, "ERROR: INVALID EVENT FOUND!");
    }
  }

  void apply_event_two_constrained_pvertices_meet(
    const PVertex& /* pvertex */,
    const PVertex& /* pother  */,
    const Event&   /* event   */) {

    CGAL_assertion_msg(false,
    "TODO: IMPLEMENT TWO CONSTRAINED PVERTICES MEET EVENT!");
  }

  void apply_event_two_unconstrained_pvertices_meet(
    const PVertex& /* pvertex */,
    const PVertex& /* pother  */,
    const Event&   /* event   */) {

    CGAL_assertion_msg(false,
    "ERROR: TWO UNCONSTRAINED PVERTICES MEET! DO WE HAVE A CONCAVE POLYGON?");
  }

  void apply_event_constrained_pvertex_meets_free_pvertex(
    const PVertex& pvertex,
    const PVertex& pother,
    const Event& event) {

    CGAL_assertion(m_data.has_iedge(pvertex));
    if (m_data.transfer_pvertex_via_iedge(pvertex, pother)) {

      if (m_data.has_iedge(pvertex)) {
        remove_events(m_data.iedge(pvertex), pvertex.first);
      }
      if (m_data.has_iedge(pother)) {
        remove_events(m_data.iedge(pother) , pother.first);
      }
      compute_events_of_pvertices(
        event.time(), std::array<PVertex, 2>{pvertex, pother});

      PVertex prev, next;
      std::tie(prev, next) = m_data.border_prev_and_next(pvertex);
      PVertex pthird = prev;
      if (pthird == pother) {
        pthird = next;
      } else {
        CGAL_assertion(next == pother);
      }

      if (m_data.has_iedge(pthird)) {
        remove_events(m_data.iedge(pthird), pthird.first);
      }
      compute_events_of_pvertices(
        event.time(), std::array<PVertex, 1>{pthird});

    } else {

      if (m_data.has_iedge(pvertex)) {
        remove_events(m_data.iedge(pvertex), pvertex.first);
      }
      compute_events_of_pvertices(
        event.time(), std::array<PVertex, 1>{pvertex});
    }
  }

  void apply_event_constrained_pvertex_meets_iedge(
    const PVertex& /* pvertex */,
    const IEdge& /* iedge */,
    const Event& /* event */) {

    CGAL_assertion_msg(false,
    "ERROR: CONSTRAINED PVERTEX MEETS IEDGE! WHAT IS WRONG?");
  }

  const bool apply_event_unconstrained_pedge_meets_iedge(
    const PVertex& pvertex,
    const IEdge& iedge,
    const Event& event) {

    bool is_event_happend = false;
    const auto pface = m_data.pface_of_pvertex(pvertex);

    std::vector<PFace> nfaces;
    m_data.non_null_pfaces_around_pvertex(pvertex, nfaces);
    CGAL_assertion(nfaces.size() == 1);
    CGAL_assertion(nfaces[0] == pface);

    const auto prev = m_data.prev(pvertex);
    const auto next = m_data.next(pvertex);
    const auto isegment = m_data.segment_2(pvertex.first, iedge);

    for (const auto& pother : { prev, next }) {
      const Segment_2 segment(
        m_data.point_2(pother , event.time()),
        m_data.point_2(pvertex, event.time()));
      CGAL_assertion(segment.squared_length() != FT(0));

      bool both_are_free = true;
      if (m_data.has_iedge(pvertex) || m_data.has_iedge(pother)) {
        both_are_free = false;
      }

      if (both_are_free && KSR::are_parallel(segment, isegment)) {
        remove_events(pvertex);
        remove_events(pother);

        m_data.non_null_pfaces_around_pvertex(pother, nfaces);
        CGAL_assertion(nfaces.size() == 1);
        CGAL_assertion(nfaces[0] == pface);

        bool collision, bbox_reached;
        std::tie(collision, bbox_reached) = m_data.collision_occured(pvertex, iedge);
        // std::tie(collision, bbox_reached) = m_data.is_occupied(pvertex, iedge);

        bool collision_other, bbox_reached_other;
        std::tie(collision_other, bbox_reached_other) = m_data.collision_occured(pother, iedge);
        // std::tie(collision_other, bbox_reached_other) = m_data.is_occupied(pother, iedge);

        if (m_debug) {
          std::cout << "- collision/bbox: " << collision << "/" << bbox_reached << std::endl;
          std::cout << "- other/bbox: " << collision_other << "/" << bbox_reached_other << std::endl;
          std::cout << "- k intersections befor: " << m_data.k(pface) << std::endl;
        }

        bool stop = false;
        if (bbox_reached) {

          CGAL_assertion(bbox_reached_other);
          if (m_debug) std::cout << "- pv po bbox" << std::endl;
          stop = true;

        } else if (bbox_reached_other) {

          CGAL_assertion(bbox_reached);
          if (m_debug) std::cout << "- po pv bbox" << std::endl;
          stop = true;

        } else if ((collision || collision_other) && m_data.k(pface) == 1) {

          if (m_debug) std::cout << "- pv po k stop" << std::endl;
          stop = true;

        } else if ((collision || collision_other) && m_data.k(pface) > 1) {

          if (m_debug) std::cout << "- pv po k continue" << std::endl;
          m_data.k(pface)--;

        } else {

          CGAL_assertion(!collision && !collision_other);
          if (m_debug) std::cout << "- pv po continue" << std::endl;
          if (
            m_data.is_occupied(pvertex, iedge).first ||
            m_data.is_occupied(pother , iedge).first) {

            CGAL_assertion_msg(false,
            "ERROR: TWO PVERTICES SNEAK TO THE OTHER SIDE EVEN WHEN WE HAVE A POLYGON!");
          }
        }

        CGAL_assertion(m_data.k(pface) >= 1);
        if (m_debug) {
          // std::cout << "PFACE: " << m_data.centroid_of_pface(pface) << std::endl;
          std::cout << "- k intersections after: " << m_data.k(pface) << std::endl;
        }

        if (stop) { // polygon stops
          m_data.crop_pedge_along_iedge(pvertex, pother, iedge);
          remove_events(iedge, pvertex.first);
          compute_events_of_pvertices(
            event.time(), std::array<PVertex, 2>{pvertex, pother});
        } else { // polygon continues beyond the edge
          PVertex pv0, pv1;
          std::tie(pv0, pv1) =
            m_data.propagate_pedge_beyond_iedge(pvertex, pother, iedge, m_data.k(pface));
          remove_events(iedge, pvertex.first);
          compute_events_of_pvertices(
            event.time(), std::array<PVertex, 4>{pvertex, pother, pv0, pv1});
        }

        CGAL_assertion(m_data.has_iedge(pvertex));
        CGAL_assertion(m_data.has_iedge(pother));
        CGAL_assertion(m_data.iedge(pvertex) == m_data.iedge(pother));
        is_event_happend = true;
        break;
      }
    }
    return is_event_happend;
  }

  void apply_event_unconstrained_pvertex_meets_iedge(
    const PVertex& pvertex,
    const IEdge& iedge,
    const Event& event) {

    CGAL_assertion(!m_data.has_iedge(pvertex));
    const auto pface = m_data.pface_of_pvertex(pvertex);
    remove_events(pvertex);

    std::vector<PFace> nfaces;
    m_data.non_null_pfaces_around_pvertex(pvertex, nfaces);
    CGAL_assertion(nfaces.size() == 1);
    CGAL_assertion(nfaces[0] == pface);

    bool collision, bbox_reached;
    std::tie(collision, bbox_reached) = m_data.collision_occured(pvertex, iedge);
    // std::tie(collision, bbox_reached) = m_data.is_occupied(pvertex, iedge);

    if (m_debug) {
      std::cout << "- collision/bbox: " << collision << "/" << bbox_reached << std::endl;
      std::cout << "- k intersections befor: " << m_data.k(pface) << std::endl;
    }

    bool stop = false;
    if (bbox_reached) {

      if (m_debug) std::cout << "- pv k bbox" << std::endl;
      stop = true;

    } else if (collision && m_data.k(pface) == 1) {

      if (m_debug) std::cout << "- pv k stop" << std::endl;
      stop = true;

    } else if (collision && m_data.k(pface) > 1) {

      if (m_debug) std::cout << "- pv k continue" << std::endl;
      m_data.k(pface)--;

    } else {
      if (m_debug) std::cout << "- pv continue" << std::endl;
    }

    CGAL_assertion(m_data.k(pface) >= 1);
    if (m_debug) {
      // std::cout << "PFACE: " << m_data.centroid_of_pface(pface) << std::endl;
      std::cout << "- k intersections after: " << m_data.k(pface) << std::endl;
    }

    if (stop) { // polygon stops
      const PVertex pother =
        m_data.crop_pvertex_along_iedge(pvertex, iedge);
      remove_events(iedge, pvertex.first);
      compute_events_of_pvertices(
        event.time(), std::array<PVertex, 2>{pvertex, pother});
    } else { // polygon continues beyond the edge
      const std::array<PVertex, 3> pvertices =
        m_data.propagate_pvertex_beyond_iedge(pvertex, iedge, m_data.k(pface));
      remove_events(iedge, pvertex.first);
      compute_events_of_pvertices(event.time(), pvertices);
    }
    CGAL_assertion(m_data.has_iedge(pvertex));
  }

  void apply_event_constrained_pvertex_meets_ivertex(
    const PVertex& pvertex,
    const IVertex& ivertex,
    const Event& event) {

    // First, let's gather all pvertices that will get merged.
    const std::vector<PVertex> crossed_pvertices =
      m_data.pvertices_around_ivertex(pvertex, ivertex);

    if (m_debug) {
      std::cout << "- found " << crossed_pvertices.size() <<
      " pvertices ready to be merged: " << std::endl;
      for (const auto& crossed_pvertex : crossed_pvertices) {
        std::cout << m_data.point_3(crossed_pvertex) << std::endl;
      }
    }

    // Remove associated events.
    for (std::size_t i = 1; i < crossed_pvertices.size() - 1; ++i) {
      remove_events(crossed_pvertices[i]);
    }

    // Merge them and get the newly created pvertices.
    CGAL_assertion(!m_data.has_ivertex(pvertex));
    std::vector<IEdge> crossed_iedges;
    const std::vector<PVertex> new_pvertices =
      m_data.merge_pvertices_on_ivertex(
        m_min_time, m_max_time, pvertex, ivertex, crossed_pvertices, crossed_iedges);

    // Remove all events of the crossed iedges.
    for (const auto& crossed_iedge : crossed_iedges) {
      remove_events(crossed_iedge, pvertex.first);
    }

    // And compute new events.
    CGAL_assertion(new_pvertices.size() > 0);
    compute_events_of_pvertices(event.time(), new_pvertices);
  }

  void apply_event_unconstrained_pvertex_meets_ivertex(
    const PVertex& /* pvertex */,
    const IVertex& /* ivertex */,
    const Event&   /* event   */) {

    CGAL_assertion_msg(false,
    "TODO: IMPLEMENT UNCONSTRAINED PVERTEX MEETS IVERTEX EVENT!");
  }

  void remove_events(const IEdge& iedge, const KSR::size_t support_plane_idx) {
    m_queue.erase_vertex_events(iedge, support_plane_idx);
    // std::cout << "erasing events for iedge: " << m_data.str(iedge) << std::endl;
    // std::cout << "iedge: " << m_data.segment_3(iedge) << std::endl;
  }

  void remove_events(const PVertex& pvertex) {
    m_queue.erase_vertex_events(pvertex);
    // std::cout << "erasing events for pvertex: " << m_data.str(pvertex) << std::endl;
    // std::cout << "pvertex: " << m_data.point_3(pvertex) << std::endl;
  }
};

} // namespace CGAL

#endif // CGAL_KINETIC_SHAPE_RECONSTRUCTION_3_H
