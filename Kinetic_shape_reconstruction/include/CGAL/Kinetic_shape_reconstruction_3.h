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
#include <CGAL/Real_timer.h>

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
  using Timer        = CGAL::Real_timer;

private:
  const bool m_debug;
  const bool m_verbose;
  Event_queue m_queue;
  FT m_min_time;
  FT m_max_time;
  Initializer m_initializer;
  Data_structure m_data;
  std::size_t m_num_events;

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
  m_data(m_debug),
  m_num_events(0)
  { }

  template<
  typename InputRange,
  typename PolygonMap,
  typename NamedParameters>
  const bool partition(
    const InputRange& input_range,
    const PolygonMap polygon_map,
    const NamedParameters& np) {

    Timer timer;
    const unsigned int k = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::k_intersections), 1);
    unsigned int n = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::n_subdivisions), 0);
    FT enlarge_bbox_ratio = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::enlarge_bbox_ratio), FT(11) / FT(10));
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

    if (enlarge_bbox_ratio < FT(1)) {
      CGAL_warning_msg(enlarge_bbox_ratio >= FT(1),
      "WARNING: YOU SET ENLARGE_BBOX_RATIO < 1.0! THE VALID RANGE IS [1.0, +INF). SETTING TO 1.0!");
      enlarge_bbox_ratio = FT(1);
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

    timer.reset();
    timer.start();
    const FT time_step = static_cast<FT>(m_initializer.initialize(
      input_range, polygon_map, k, CGAL::to_double(enlarge_bbox_ratio), reorient));
    m_initializer.convert(m_data);
    m_data.set_limit_lines();
    m_data.precompute_iedge_data();
    CGAL_assertion(m_data.check_integrity());
    timer.stop();
    const double time_to_initialize = timer.time();

    if (k == 0) {
      CGAL_warning_msg(k > 0,
      "WARNING: YOU SET K TO 0! THAT MEANS NO PROPAGATION! THE VALID VALUES ARE {1,2,...}. INTERSECT AND RETURN!");
      return false;
    }

    // if (m_verbose) {
    //   std::cout << std::endl << "* initialization (sec.): " << time_to_initialize << std::endl;
    //   std::cout << "POLYGON SPLITTER SUCCESS!" << std::endl << std::endl;
    // }
    // exit(EXIT_SUCCESS);

    // Output planes.
    // for (std::size_t i = 6; i < m_data.number_of_support_planes(); ++i) {
    //   std::cout << m_data.support_plane(i).plane() << std::endl;
    // }

    if (m_verbose) {
      std::cout << std::endl << "--- RUNNING THE QUEUE:" << std::endl;
      std::cout << "* propagation started" << std::endl;
    }

    timer.reset();
    timer.start();
    std::size_t num_iterations = 0;
    m_min_time = FT(0);
    m_max_time = time_step;
    CGAL_assertion(m_min_time >= FT(0) && m_max_time >= m_min_time);
    std::size_t global_iteration = 0;
    while (initialize_queue()) {

      global_iteration = run(k, global_iteration);
      m_min_time = m_max_time;
      m_max_time += time_step;
      CGAL_assertion(m_data.check_integrity());
      ++num_iterations;

      if (m_verbose && !m_debug) {
        if ((num_iterations % 50) == 0) {
          std::cout << ".................................................." << std::endl;
        }
      }

      if (num_iterations > 100000000) {
        CGAL_assertion_msg(false, "DEBUG ERROR: WHY SO MANY ITERATIONS?");
        return false;
      }
    }
    timer.stop();
    const double time_to_propagate = timer.time();

    if (m_verbose) {
      std::cout << "* propagation finished" << std::endl;
      std::cout << "* number of iterations: " << num_iterations   << std::endl;
      std::cout << "* number of events: "     << global_iteration << std::endl;
    }
    m_num_events = global_iteration;

    timer.reset();
    timer.start();
    if (m_verbose) std::cout << std::endl << "--- FINALIZING PARTITION:" << std::endl;
    if (m_debug) dump(m_data, "jiter-final-a-result");
    m_data.finalize();
    if (m_verbose) std::cout << "* checking final mesh integrity ...";
    CGAL_assertion(m_data.check_integrity(true, true, true));
    if (m_verbose) std::cout << " done" << std::endl;
    if (m_debug) dump(m_data, "jiter-final-b-result");
    // std::cout << std::endl << "CLEANING SUCCESS!" << std::endl << std::endl;
    // exit(EXIT_SUCCESS);
    if (m_verbose) std::cout << "* getting volumes ..." << std::endl;
    m_data.create_polyhedra();
    timer.stop();
    const double time_to_finalize = timer.time();
    if (m_verbose) {
      std::cout << "* found " << m_data.number_of_volumes(-1) << " volumes" << std::endl;
    }

    if (m_verbose) std::cout << std::endl << "--- TIMING (sec.):" << std::endl;
    const double total_time =
      time_to_initialize + time_to_propagate + time_to_finalize;
    if (m_verbose) {
      std::cout << "* initialization: " << time_to_initialize << std::endl;
      std::cout << "* propagation: "    << time_to_propagate  << std::endl;
      std::cout << "* finalization: "   << time_to_finalize   << std::endl;
      std::cout << "* total time: "     << total_time         << std::endl;
    }
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
      input_range, point_map, normal_map, semantic_map, m_data, m_verbose, m_debug);
    bool success = reconstruction.detect_planar_shapes(np);
    if (!success) {
      CGAL_assertion_msg(false, "ERROR: RECONSTRUCTION, DETECTING PLANAR SHAPES FAILED!");
      return false;
    }
    // exit(EXIT_SUCCESS);

    success = reconstruction.regularize_planar_shapes(np);
    if (!success) {
      CGAL_assertion_msg(false, "ERROR: RECONSTRUCTION, REGULARIZATION FAILED!");
      return false;
    }
    // exit(EXIT_SUCCESS);

    success = partition(
      reconstruction.planar_shapes(), reconstruction.polygon_map(), np);
    if (!success) {
      CGAL_assertion_msg(false, "ERROR: RECONSTRUCTION, PARTITION FAILED!");
      return false;
    }
    // exit(EXIT_SUCCESS);

    success = reconstruction.compute_model(np);
    if (!success) {
      CGAL_assertion_msg(false, "ERROR: RECONSTRUCTION, COMPUTING MODEL FAILED!");
      return false;
    }
    return success;
  }

  /*******************************
  **         STATISTICS         **
  ********************************/

  const std::size_t number_of_events() const {
    return m_num_events;
  }

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
    const std::size_t sp_idx = static_cast<std::size_t>(support_plane_idx);
    return static_cast<std::size_t>(m_data.mesh(sp_idx).number_of_vertices());
  }

  const std::size_t number_of_edges(const int support_plane_idx = -1) const {

    CGAL_assertion(support_plane_idx < number_of_support_planes());
    if (support_plane_idx >= number_of_support_planes()) return std::size_t(-1);
    if (support_plane_idx < 0) {
      return m_data.igraph().number_of_edges();
    }

    CGAL_assertion(support_plane_idx >= 0);
    const std::size_t sp_idx = static_cast<std::size_t>(support_plane_idx);
    return static_cast<std::size_t>(m_data.mesh(sp_idx).number_of_edges());
  }

  const std::size_t number_of_faces(const int support_plane_idx = -1) const {

    CGAL_assertion(support_plane_idx < number_of_support_planes());
    if (support_plane_idx >= number_of_support_planes()) return std::size_t(-1);
    if (support_plane_idx < 0) {
      std::size_t num_all_faces = 0;
      for (int i = 0; i < number_of_support_planes(); ++i) {
        const std::size_t num_faces = static_cast<std::size_t>(
          m_data.mesh(static_cast<std::size_t>(i)).number_of_faces());
        num_all_faces += num_faces;
      }
      return num_all_faces;
    }

    CGAL_assertion(support_plane_idx >= 0);
    const std::size_t sp_idx = static_cast<std::size_t>(support_plane_idx);
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
    const std::size_t sp_idx = static_cast<std::size_t>(support_plane_idx);
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
    const std::size_t sp_idx = static_cast<std::size_t>(support_plane_idx);
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
    FaceOutputIterator faces,
    const int support_plane_idx = -1,
    const int begin = 0) const {

    KSR::Indexer<IVertex> indexer;
    CGAL_assertion(support_plane_idx < number_of_support_planes());
    if (support_plane_idx >= number_of_support_planes()) return faces;
    if (support_plane_idx < 0) {
      const auto all_ivertices = m_data.ivertices();
      for (const auto ivertex : all_ivertices) indexer(ivertex);
      for (int i = begin; i < number_of_support_planes(); ++i) {
        const std::size_t sp_idx = static_cast<std::size_t>(i);
        output_partition_faces(faces, indexer, sp_idx);
      }
      return faces;
    }

    CGAL_assertion(support_plane_idx >= 0);
    const std::size_t sp_idx = static_cast<std::size_t>(support_plane_idx);
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
    const std::size_t sp_idx = static_cast<std::size_t>(support_plane_idx);

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

    const auto& model = m_data.reconstructed_model();
    CGAL_assertion(model.pfaces.size() > 0);

    std::size_t num_vertices = 0;
    KSR::Indexer<IVertex> indexer;

    std::vector<std::size_t> face;
    const auto& pfaces = model.pfaces;
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

  void output_reconstructed_model(Polygon_mesh& polygon_mesh) const {

    std::vector<Point_3> vertices;
    std::vector< std::vector<std::size_t> > faces;
    output_reconstructed_model(
      std::back_inserter(vertices), std::back_inserter(faces));
    CGAL::Polygon_mesh_processing::orient_polygon_soup(vertices, faces);
    polygon_mesh.clear();
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(
      vertices, faces, polygon_mesh);
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
    m_num_events = 0;
  }

private:

  template<typename FaceOutputIterator>
  FaceOutputIterator output_partition_faces(
    FaceOutputIterator faces, KSR::Indexer<IVertex>& indexer,
    const std::size_t sp_idx) const {

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
    for (std::size_t i = 0; i < m_data.number_of_support_planes(); ++i) {
      const auto& iedges   = m_data.iedges(i);
      const auto& segments = m_data.isegments(i);
      const auto& bboxes   = m_data.ibboxes(i);
      for (const auto pvertex : m_data.pvertices(i)) {
        if (compute_events_of_pvertex(pvertex, iedges, segments, bboxes)) {
          still_running = true;
        }
      }
    }
    m_data.update_positions(m_min_time);
    return still_running;
  }

  const bool compute_events_of_pvertex(
    const PVertex& pvertex,
    const std::vector<IEdge>& iedges,
    const std::vector<Segment_2>& segments,
    const std::vector<Bbox_2>& bboxes) {

    CGAL_assertion(iedges.size() > 0);
    CGAL_assertion(iedges.size() == segments.size());
    CGAL_assertion(iedges.size() == bboxes.size());

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
    const PVertex& pvertex, const Segment_2& pv_segment, const Bbox_2& pv_bbox) {

    // const bool is_event_found =
    // try_pvertices_to_ivertex_event(pvertex, pv_segment, pv_bbox);
    // if (!is_event_found) return;

    try_pvertex_to_pvertex_constrained_event(pvertex, pv_segment, pv_bbox);
    try_pvertex_to_ivertex_constrained_event(pvertex, pv_segment, pv_bbox);
  }

  const bool try_pvertices_to_ivertex_event(
    const PVertex& pvertex, const Segment_2& pv_segment, const Bbox_2& pv_bbox) {
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

      CGAL_assertion(m_data.has_iedge(pvertex));
      const auto iedge = m_data.iedge(pvertex);

      const auto isource = m_data.source(iedge);
      const auto itarget = m_data.target(iedge);

      const auto source = m_data.point_2(pvertex.first, isource);
      const auto target = m_data.point_2(pvertex.first, itarget);

      const FT tol = KSR::tolerance<FT>();
      const FT dist1 = KSR::distance(inter, source);
      const FT dist2 = KSR::distance(inter, target);

      // std::cout << "tol: " << tol << std::endl;
      // std::cout << "dist 1: " << dist1 << std::endl;
      // std::cout << "dist 2: " << dist2 << std::endl;

      Point_2 ipoint;
      IVertex ivertex = m_data.null_ivertex();
      if (dist1 < tol) {
        CGAL_assertion(dist2 >= tol);
        ipoint = source; ivertex = isource;
      } else if (dist2 < tol) {
        CGAL_assertion(dist1 >= tol);
        ipoint = target; ivertex = itarget;
      }

      if (ivertex != m_data.null_ivertex()) {
        CGAL_assertion(ipoint != Point_2());

        const auto& pinit = pv_segment.source();
        const FT distance = KSR::distance(pinit, ipoint);
        const FT time = distance / m_data.speed(pvertex);

        // Should I break here?
        is_event_found = true;
        CGAL_assertion(time < m_max_time - m_min_time);
        m_queue.push(Event(true, pvertex, pother, ivertex, m_min_time + time));
        CGAL_assertion_msg(false, "TODO: TRY PVERTICES TO IVERTEX EVENT!");
      }
    }
    return is_event_found;
  }

  void try_pvertex_to_pvertex_constrained_event(
    const PVertex& pvertex, const Segment_2& pv_segment, const Bbox_2& pv_bbox) {

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

      const auto& pinit = pv_segment.source();
      const FT distance = KSR::distance(pinit, inter);
      const FT time = distance / m_data.speed(pvertex);

      // Constrained pvertex to another pvertex event.
      CGAL_assertion(time < m_max_time - m_min_time);
      m_queue.push(Event(true, pvertex, pother, m_min_time + time));

      // std::cout << "pvertex: " << m_data.point_3(pvertex) << std::endl;
      // std::cout << "pother: "  << m_data.point_3(pother)  << std::endl;
    }
  }

  void try_pvertex_to_ivertex_constrained_event(
    const PVertex& pvertex, const Segment_2& pv_segment, const Bbox_2& pv_bbox) {

    CGAL_assertion(m_data.has_iedge(pvertex));
    const auto iedge = m_data.iedge(pvertex);
    for (const auto& ivertex : { m_data.source(iedge), m_data.target(iedge) }) {
      if (!m_data.is_active(ivertex)) {
        continue;
      }

      const Point_2 ipoint = m_data.to_2d(pvertex.first, ivertex);
      const auto vec1 = pv_segment.to_vector();
      const auto& pinit = pv_segment.source();
      const Vector_2 vec2(pinit, ipoint);
      const FT dot_product = vec1 * vec2;
      if (dot_product < FT(0)) { // opposite directions
        continue;
      }

      const FT distance = KSR::distance(pinit, ipoint);
      const FT time = distance / m_data.speed(pvertex);

      // Constrained pvertex to ivertex event.
      if (time < m_max_time - m_min_time) {

        CGAL_assertion(time < m_max_time - m_min_time);
        m_queue.push(Event(true, pvertex, ivertex, m_min_time + time));

        // std::cout << "pvertex: " << m_data.point_3(pvertex) << std::endl;
        // std::cout << "ivertex: " << m_data.point_3(ivertex) << std::endl;
      }
    }
  }

  void compute_events_of_unconstrained_pvertex(
    const PVertex& pvertex,
    const Segment_2& pv_segment,
    const Bbox_2& pv_bbox,
    const std::vector<IEdge>& iedges,
    const std::vector<Segment_2>& segments,
    const std::vector<Bbox_2>& bboxes) {

    try_pvertex_to_iedge_unconstrained_event(
      pvertex, pv_segment, pv_bbox, iedges, segments, bboxes);
  }

  void try_pvertex_to_iedge_unconstrained_event(
    const PVertex& pvertex,
    const Segment_2& pv_segment,
    const Bbox_2& pv_bbox,
    const std::vector<IEdge>& iedges,
    const std::vector<Segment_2>& segments,
    const std::vector<Bbox_2>& bboxes) {

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
      const auto& pinit = pv_segment.source();
      // const bool is_event_found = try_pvertex_to_ivertex_unconstrained_event(
      //   pvertex, iedge, inter, pinit);

      // Otherwise we add unconstrained pvertex to iedge event.
      // if (!is_event_found) {

      const FT distance = KSR::distance(pinit, inter);
      const FT time = distance / m_data.speed(pvertex);
      CGAL_assertion(time < m_max_time - m_min_time);
      m_queue.push(Event(false, pvertex, iedge, m_min_time + time));

      // }

      // std::cout << "pvertex: " << m_data.point_3(pvertex) << std::endl;
      // std::cout << "iedge: "   << m_data.segment_3(iedge) << std::endl;
    }
  }

  const bool try_pvertex_to_ivertex_unconstrained_event(
    const PVertex& pvertex, const IEdge& iedge,
    const Point_2& inter, const Point_2& pinit) {

    bool is_event_found = false;
    const auto isource = m_data.source(iedge);
    const auto itarget = m_data.target(iedge);

    const auto source = m_data.point_2(pvertex.first, isource);
    const auto target = m_data.point_2(pvertex.first, itarget);

    const FT tol = KSR::tolerance<FT>();
    const FT dist1 = KSR::distance(inter, source);
    const FT dist2 = KSR::distance(inter, target);

    // std::cout << "tol: " << tol << std::endl;
    // std::cout << "dist 1: " << dist1 << std::endl;
    // std::cout << "dist 2: " << dist2 << std::endl;

    Point_2 ipoint;
    IVertex ivertex = m_data.null_ivertex();
    if (dist1 < tol) {
      CGAL_assertion(dist2 >= tol);
      ipoint = source; ivertex = isource;
    } else if (dist2 < tol) {
      CGAL_assertion(dist1 >= tol);
      ipoint = target; ivertex = itarget;
    }

    if (ivertex != m_data.null_ivertex()) {
      CGAL_assertion(ipoint != Point_2());
      const FT distance = KSR::distance(pinit, ipoint);
      const FT time = distance / m_data.speed(pvertex);
      CGAL_assertion(time < m_max_time - m_min_time);
      m_queue.push(Event(false, pvertex, ivertex, m_min_time + time));
      is_event_found = true;
    }

    // CGAL_assertion_msg(false, "TODO: ADD PVERTEX TO IVERTEX UNCONSTRAINED EVENT!");
    return is_event_found;
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
        // const std::size_t sp_debug_idx = 23;
        // dump_2d_surface_mesh(m_data, sp_debug_idx, "iter-" + std::to_string(iteration) +
        // "-surface-mesh-" + std::to_string(sp_debug_idx));
      }

      m_data.update_positions(current_time);
      if (m_debug) {
        std::cout << std::endl << "* APPLYING " << iteration << ": " << event << std::endl;
      }
      ++iteration;

      // if (iteration == 80) {
      //   exit(EXIT_FAILURE);
      // }

      apply(event, k);
      CGAL_assertion(m_data.check_integrity());
    }
    return iteration;
  }

  void apply(const Event& event, const unsigned int k) {

    const auto pvertex = event.pvertex();
    if (event.is_pvertices_to_ivertex()) {

      const auto pother  = event.pother();
      const auto ivertex = event.ivertex();
      apply_event_pvertices_meet_ivertex(pvertex, pother, ivertex, event);

    } else if (event.is_pvertex_to_pvertex()) {
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

  // INVALID EVENTS!
  void apply_event_two_unconstrained_pvertices_meet(
    const PVertex& /* pvertex */,
    const PVertex& /* pother  */,
    const Event&   /* event   */) {

    CGAL_assertion_msg(false,
    "ERROR: TWO UNCONSTRAINED PVERTICES MEET! DO WE HAVE A CONCAVE POLYGON?");
  }

  void apply_event_two_constrained_pvertices_meet(
    const PVertex& /* pvertex */,
    const PVertex& /* pother  */,
    const Event&   /* event   */) {

    CGAL_assertion_msg(false,
    "ERROR: TWO CONSTRAINED PVERTICES MEET! CAN IT HAPPEN?");
  }

  void apply_event_constrained_pvertex_meets_iedge(
    const PVertex& /* pvertex */,
    const IEdge&   /* iedge   */,
    const Event&   /* event   */) {

    CGAL_assertion_msg(false,
    "ERROR: CONSTRAINED PVERTEX MEETS IEDGE! WHAT IS WRONG?");
  }

  void apply_event_pvertices_meet_ivertex(
    const PVertex& pvertex, const PVertex& pother,
    const IVertex& /* ivertex */, const Event& /* event */) {

    CGAL_assertion( m_data.has_iedge(pvertex));
    CGAL_assertion(!m_data.has_iedge(pother));
    CGAL_assertion_msg(false,
    "ERROR: PVERTICES MEET IVERTEX! IT SHOULD NOT EVER HAPPEN!");
  }

  void apply_event_unconstrained_pvertex_meets_ivertex(
    const PVertex& pvertex, const IVertex& ivertex, const Event& event) {

    CGAL_assertion(!m_data.has_iedge(pvertex));
    CGAL_assertion( m_data.has_one_pface(pvertex));

    CGAL_assertion_msg(false,
    "ERROR: UNCONSTRAINED PVERTEX MEETS IVERTEX! IT SHOULD NOT EVER HAPPEN!");
    apply_event_pvertex_meets_ivertex(pvertex, ivertex, event);
  }

  // VALID EVENTS!
  void apply_event_pvertex_meets_ivertex(
    const PVertex& pvertex, const IVertex& ivertex, const Event& event) {

    // First, let's gather all pvertices that will get merged.
    const std::vector<PVertex> crossed_pvertices =
      m_data.pvertices_around_ivertex(pvertex, ivertex);

    // Remove associated events.
    CGAL_assertion(crossed_pvertices.size() >= 3);
    for (std::size_t i = 1; i < crossed_pvertices.size() - 1; ++i) {
      remove_events(crossed_pvertices[i]);
    }

    // Merge them and get the newly created pvertices.
    CGAL_assertion(!m_data.has_ivertex(pvertex));
    std::vector< std::pair<IEdge, bool> > crossed_iedges;
    const std::vector<PVertex> pvertices =
      m_data.merge_pvertices_on_ivertex(
        m_min_time, m_max_time, ivertex, crossed_pvertices, crossed_iedges);

    // Remove all events of the crossed iedges.
    CGAL_assertion(crossed_iedges.size() >= 1);
    for (const auto& crossed_iedge : crossed_iedges) {
      // TODO: SHOULD I LEAVE THIS CHECK? WILL IT MAKE THE CODE FASTER?
      // if (crossed_iedges[ip].second) {
        // bla bla
      // }
      const auto& iedge = crossed_iedge.first;
      remove_events(iedge, pvertex.first);
    }

    // And compute new events.
    CGAL_assertion(pvertices.size() > 0);
    compute_events_of_pvertices(event.time(), pvertices);
    // CGAL_assertion_msg(false, "TODO: PVERTEX MEETS IVERTEX!");
  }

  void apply_event_unconstrained_pvertex_meets_iedge(
    const PVertex& pvertex, const IEdge& iedge, const Event& event) {

    CGAL_assertion(!m_data.has_iedge(pvertex));
    CGAL_assertion( m_data.has_one_pface(pvertex));

    remove_events(pvertex);
    const bool stop = check_stop_condition(pvertex, iedge);

    if (stop) { // polygon stops
      const PVertex pother =
        m_data.crop_pvertex_along_iedge(pvertex, iedge);
      const std::array<PVertex, 2> pvertices = {pvertex, pother};
      remove_events(iedge, pvertex.first);
      compute_events_of_pvertices(event.time(), pvertices);
    } else { // polygon continues beyond the iedge
      const std::array<PVertex, 3> pvertices =
        m_data.propagate_pvertex_beyond_iedge(pvertex, iedge);
      remove_events(iedge, pvertex.first);
      compute_events_of_pvertices(event.time(), pvertices);
    }
    CGAL_assertion(m_data.has_iedge(pvertex));
    // CGAL_assertion_msg(false, "TODO: UNCONSTRAINED PVERTEX MEETS IEDGE!");
  }

  const bool apply_event_unconstrained_pedge_meets_iedge(
    const PVertex& pvertex, const IEdge& iedge, const Event& event) {

    bool is_event_happend = false;
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
        CGAL_assertion(!m_data.has_iedge(pother));
        CGAL_assertion(!m_data.has_iedge(pvertex));

        CGAL_assertion(m_data.has_one_pface(pother));
        CGAL_assertion(m_data.has_one_pface(pvertex));

        remove_events(pother);
        remove_events(pvertex);

        const bool stop = check_stop_condition(pvertex, pother, iedge);

        if (stop) { // polygon stops
          m_data.crop_pedge_along_iedge(pvertex, pother, iedge);
          const auto pvertices = std::array<PVertex, 2>{pvertex, pother};
          remove_events(iedge, pvertex.first);
          compute_events_of_pvertices(event.time(), pvertices);
        } else { // polygon continues beyond the edge
          PVertex pv0, pv1;
          std::tie(pv0, pv1) = m_data.propagate_pedge_beyond_iedge(pvertex, pother, iedge);
          const auto pvertices = std::array<PVertex, 4>{pvertex, pother, pv0, pv1};
          remove_events(iedge, pvertex.first);
          compute_events_of_pvertices(event.time(), pvertices);
        }

        CGAL_assertion(m_data.has_iedge(pother));
        CGAL_assertion(m_data.has_iedge(pvertex));
        CGAL_assertion(m_data.iedge(pvertex) == m_data.iedge(pother));
        is_event_happend = true;
        // CGAL_assertion_msg(false, "TODO: UNCONSTRAINED PEDGE MEETS IEDGE!");
        break;
      }
    }
    return is_event_happend;
  }

  void apply_event_constrained_pvertex_meets_ivertex(
    const PVertex& pvertex, const IVertex& ivertex, const Event& event) {

    CGAL_assertion(m_data.has_iedge(pvertex));
    apply_event_pvertex_meets_ivertex(pvertex, ivertex, event);
    // CGAL_assertion_msg(false, "TODO: CONSTRAINED PVERTEX MEETS IVERTEX!");
  }

  void apply_event_constrained_pvertex_meets_free_pvertex(
    const PVertex& pvertex, const PVertex& pother, const Event& event) {

    CGAL_assertion( m_data.has_iedge(pvertex));
    CGAL_assertion(!m_data.has_iedge(pother));

    if (m_data.transfer_pvertex_via_iedge(pvertex, pother)) {

      // Check the first two pvertices.
      if (m_data.has_iedge(pvertex)) {
        remove_events(m_data.iedge(pvertex), pvertex.first);
      }
      if (m_data.has_iedge(pother)) {
        remove_events(m_data.iedge(pother) , pother.first);
      }
      const auto pvertices1 = std::array<PVertex, 2>{pvertex, pother};
      compute_events_of_pvertices(event.time(), pvertices1);

      // Check the last pvertex.
      PVertex prev, next;
      std::tie(prev, next) = m_data.border_prev_and_next(pvertex);
      PVertex pthird = prev;
      if (pthird == pother) {
        pthird = next;
      } else { CGAL_assertion(pother == next); }

      if (m_data.has_iedge(pthird)) {
        remove_events(m_data.iedge(pthird), pthird.first);
      }
      const auto pvertices2 = std::array<PVertex, 1>{pthird};
      compute_events_of_pvertices(event.time(), pvertices2);

    } else {

      if (m_data.has_iedge(pvertex)) {
        remove_events(m_data.iedge(pvertex), pvertex.first);
      }
      const auto pvertices = std::array<PVertex, 1>{pvertex};
      compute_events_of_pvertices(event.time(), pvertices);
    }
    // CGAL_assertion_msg(false, "TODO: CONSTRAINED PVERTEX MEETS FREE PVERTEX!");
  }

  // STOP CONDITIONS!
  const bool check_stop_condition(
    const PVertex& pvertex, const IEdge& iedge) {
    return check_pvertex_meets_iedge_global_k(pvertex, iedge);
  }

  // GLOBAL STOP CONDITIONS!
  const bool check_pvertex_meets_iedge_global_k(
    const PVertex& pvertex, const IEdge& iedge) {

    if (m_debug) {
      std::cout << "- k intersections before: " << m_data.k(pvertex.first) << std::endl;
    }

    bool is_occupied_iedge, is_bbox_reached;
    std::tie(is_occupied_iedge, is_bbox_reached) = m_data.collision_occured(pvertex, iedge);
    const bool is_limit_line = m_data.update_limit_lines_and_k(pvertex, iedge, is_occupied_iedge);

    if (m_debug) {
      std::cout << "- bbox: " << is_bbox_reached  << "; " <<
      " limit: "    << is_limit_line << "; " <<
      " occupied: " << is_occupied_iedge << std::endl;
    }

    bool stop = false;
    if (is_bbox_reached) {
      if (m_debug) std::cout << "- bbox, stop" << std::endl;
      stop = true;
    } else if (is_limit_line) {
      if (m_debug) std::cout << "- limit, stop" << std::endl;
      stop = true;
    } else {
      if (m_debug) std::cout << "- free, any k, continue" << std::endl;
      stop = false;
    }
    CGAL_assertion(m_data.k(pvertex.first) >= 1);
    if (m_debug) {
      std::cout << "- k intersections after: " << m_data.k(pvertex.first) << std::endl;
    }
    // CGAL_assertion_msg(false, "TODO: CHECK PVERTEX MEETS IVERTEX GLOBAL!");
    return stop;
  }

  const bool check_stop_condition(
    const PVertex& pvertex, const PVertex& pother, const IEdge& iedge) {
    return check_pedge_meets_iedge_global_k(pvertex, pother, iedge);
  }

  const bool check_pedge_meets_iedge_global_k(
    const PVertex& pvertex, const PVertex& pother, const IEdge& iedge) {

    if (m_debug) {
      std::cout << "- k intersections before: " << m_data.k(pvertex.first) << std::endl;
    }

    bool is_occupied_iedge_1, is_bbox_reached_1;
    std::tie(is_occupied_iedge_1, is_bbox_reached_1) = m_data.collision_occured(pvertex, iedge);
    bool is_occupied_iedge_2, is_bbox_reached_2;
    std::tie(is_occupied_iedge_2, is_bbox_reached_2) = m_data.collision_occured(pother, iedge);

    const bool is_limit_line_1 = m_data.update_limit_lines_and_k(pvertex, iedge, is_occupied_iedge_1);
    const bool is_limit_line_2 = m_data.update_limit_lines_and_k(pother , iedge, is_occupied_iedge_2);

    if (m_debug) {
      std::cout << "- bbox1: " << is_bbox_reached_1  << "; " <<
      " limit1: "    << is_limit_line_1 << "; " <<
      " occupied1: " << is_occupied_iedge_1 << std::endl;
      std::cout << "- bbox2: " << is_bbox_reached_2  << "; " <<
      " limit2: "    << is_limit_line_2 << "; " <<
      " occupied2: " << is_occupied_iedge_2 << std::endl;
    }
    CGAL_assertion(is_limit_line_1 == is_limit_line_2);
    CGAL_assertion(is_bbox_reached_1 == is_bbox_reached_2);

    bool stop = false;
    if (is_bbox_reached_1 || is_bbox_reached_2) {
      if (m_debug) std::cout << "- bbox, stop" << std::endl;
      stop = true;
    } else if (is_limit_line_1 || is_limit_line_2) {
      if (m_debug) std::cout << "- limit, stop" << std::endl;
      stop = true;
    } else {
      if (m_debug) std::cout << "- free, any k, continue" << std::endl;
      CGAL_assertion(!m_data.is_sneaking_pedge(pvertex, pother, iedge));
      stop = false;
    }

    CGAL_assertion(m_data.k(pvertex.first) >= 1);
    if (m_debug) {
      std::cout << "- k intersections after: " << m_data.k(pvertex.first) << std::endl;
    }
    // CGAL_assertion_msg(false, "TODO: CHECK PEDGE MEETS IEDGE GLOBAL!");
    return stop;
  }

  // RECOMPUTE EVENTS!
  template<typename PVertexRange>
  void compute_events_of_pvertices(
    const FT last_event_time, const PVertexRange& pvertices) {

    m_min_time = m_data.current_time();
    m_data.update_positions(m_max_time);

    const auto& pfront = pvertices.front();
    CGAL_assertion(pfront != Data_structure::null_pvertex());
    const auto& iedges   = m_data.iedges(pfront.first);
    const auto& segments = m_data.isegments(pfront.first);
    const auto& bboxes   = m_data.ibboxes(pfront.first);

    for (const auto& pvertex : pvertices) {
      if (pvertex == Data_structure::null_pvertex()) continue;
      m_data.deactivate(pvertex);
    }
    for (const auto& pvertex : pvertices) {
      if (pvertex == Data_structure::null_pvertex()) continue;
      m_data.set_last_event_time(pvertex, last_event_time);
      compute_events_of_pvertex(pvertex, iedges, segments, bboxes);
    }
    for (const auto& pvertex : pvertices) {
      if (pvertex == Data_structure::null_pvertex()) continue;
      m_data.activate(pvertex);
    }
    m_data.update_positions(m_min_time);
  }

  // REMOVE EVENTS!
  // Remove events associated with the given iedge.
  void remove_events(const IEdge& iedge, const std::size_t support_plane_idx) {
    CGAL_assertion(iedge != Data_structure::null_iedge());
    m_queue.erase_vertex_events(iedge, support_plane_idx);
    // std::cout << "erasing events for iedge: " << m_data.str(iedge) << std::endl;
    // std::cout << "iedge: " << m_data.segment_3(iedge) << std::endl;
  }

  // Remove events associated with the given pvertex.
  void remove_events(const PVertex& pvertex) {
    CGAL_assertion(pvertex != Data_structure::null_pvertex());
    m_queue.erase_vertex_events(pvertex);
    // std::cout << "erasing events for pvertex: " << m_data.str(pvertex) << std::endl;
    // std::cout << "pvertex: " << m_data.point_3(pvertex) << std::endl;
  }
};

} // namespace CGAL

#endif // CGAL_KINETIC_SHAPE_RECONSTRUCTION_3_H
