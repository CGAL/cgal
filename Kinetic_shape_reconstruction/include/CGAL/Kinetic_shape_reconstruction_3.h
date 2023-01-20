// Copyright (c) 2019 GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Simon Giraudot, Dmitry Anisimov

#ifndef CGAL_KINETIC_SHAPE_RECONSTRUCTION_3_H
#define CGAL_KINETIC_SHAPE_RECONSTRUCTION_3_H

// #include <CGAL/license/Kinetic_shape_reconstruction.h>

// Boost includes.
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Named_function_parameters.h>

// CGAL includes.
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Real_timer.h>

#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_incremental_builder.h>
#include <CGAL/draw_linear_cell_complex.h>

// Internal includes.
#include <CGAL/KSR/utils.h>
#include <CGAL/KSR/debug.h>
#include <CGAL/KSR/parameters.h>

#include <CGAL/KSR_3/Data_structure.h>
#include <CGAL/KSR_3/Reconstruction.h>
#include <CGAL/KSR_3/Initializer.h>
#include <CGAL/KSR_3/FacePropagation.h>
#include <CGAL/KSR_3/Finalizer.h>

namespace CGAL {

#ifdef DOXYGEN_NS
/*!
* \ingroup PkgKineticPartition
  \brief Creates the kinetic partitioning of the bounding box.

  \tparam Kernel
    must be a model of `Kernel`. Is used for non-critical calculations.

  \tparam IntersectionKernel
    must be a model of `Kernel`. Is used for the creation of the intersection graph. An exact kernel is suggested.
*/
template<Kernel, Intersection_Kernel = CGAL::Exact_predicates_exact_constructions_kernel>
class Kinetic_partitioning_3 {
public:
  /// \name Initialization
  /// @{
  /*!
  \brief Initializes the kinetic partitioning of the bounding box.

  \tparam InputRange
  must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

  \tparam PolygonMap
  contains index ranges to form polygons from InputRange

  \tparam NamedParameters
  a sequence of \ref bgl_namedparameters "Named Parameters"

  \param input_range
   an instance of `InputRange` with 3D points and corresponding 3D normal vectors

  \param polygon_map
  a range of polygons defined by a range of indices into input_range

  \param np
  a sequence of \ref bgl_namedparameters "Named Parameters"
  among the ones listed below

  @return
  success

  \pre `input_range.size() > 0 and polygon_map.size() > 0`

  \cgalNamedParamsBegin
    \cgalParamNBegin{reorient_bbox}
      \cgalParamDescription{Use the oriented bounding box instead of the axis-aligned bounding box.}
      \cgalParamType{bool}
      \cgalParamDefault{false}
    \cgalParamNEnd
    \cgalParamNBegin{bbox_extension}
      \cgalParamDescription{Factor for extension of the bounding box of the input data to be used for the partitioning.}
      \cgalParamType{FT}
      \cgalParamDefault{1.1}
    \cgalParamNEnd
    \cgalParamNBegin{theta}
      \cgalParamDescription{The tolerance angle to snap the planes of two input polygons into one plane.}
      \cgalParamType{FT}
      \cgalParamDefault{5}
    \cgalParamNEnd
    \cgalParamNBegin{epsilon}
      \cgalParamDescription{The tolerance distance to snap the planes of two input polygons into one plane.}
      \cgalParamType{FT}
      \cgalParamDefault{5}
    \cgalParamNEnd
  \cgalNamedParamsEnd
  */
  template<
    typename InputRange,
    typename PolygonMap,
    typename NamedParameters>
  bool initialization(
    const InputRange input_range,
    const PolygonMap polygon_map,
    const NamedParameters& np);

  /*!
  \brief Propagates the kinetic polygons in the initialized partition.

    \param k
    maximum number of allowed intersections for each input polygon before its expansion stops.

    @return
    success of kinetic partitioning.

    \pre `successful initialization`
  */
  /// @}

  bool partition(std::size_t k);

  /// \name Access
  /// @{

  /*!
  \brief Number of vertices in the kinetic partitioning.

    @return
    number of vertices.

    \pre `successful partitioning`
  */
  std::size_t number_of_vertices() const;

  /*!
  \brief Number of convex faces in the kinetic partitioning.

    @return
    number of convex faces.

    \pre `successful partitioning`
  */
  std::size_t number_of_faces() const;

  /*!
  \brief Number of convex volumes created by the kinetic partitioning.

    @return
    number of convex volumes.

    \pre `successful partitioning`
  */
  std::size_t number_of_volumes() const;

  /*!
  \brief Point vector for mapping vertex indices to positions.

    @return
    vector of points.

    \pre `successful partitioning`
  */
  const std::vector<Point_3>& vertices() const;

  /*!
  \brief Vertex indices of convex face.

    \param face_index
    index of the query face.

    @return
    vector of vertex indices.

    \pre `successful partitioning`
  */
  const std::vector<std::size_t>& vertices(std::size_t face_index) const;

  /*!
  \brief Face indices of the convex volume.

    \param volume_index
    index of the query volume.

    @return
    vector of face indices.

    \pre `successful partitioning`
  */
  const std::vector<std::size_t>& face(std::size_t volume_index) const;

  /*!
  \brief Indices of adjacent volumes. Negative indices correspond to the empty spaces behind the sides of the bounding box.

    \param face_index
    index of the query face.

    @return
    pair of adjacent volumes.

    -1 zmin
    -2 ymin
    -3 xmax
    -4 ymax
    -5 xmin
    -6 zmax

    \pre `successful partitioning`
  */
  const std::pair<int, int>& neighbors(std::size_t face_index) const;

  /*!
  \brief Retrieves the input polygon this face originates from.

    \param face_index
    index of the query face.

    @return
    index into polygon_map provided on initialization.

    \pre `successful partitioning`
  */
  const std::size_t input_polygon(std::size_t face_index) const;

  /*!
   \brief Creates a linear cell complex from the kinetic partitioning.

    \tparam LCC
    Linear_cell_complex_for_combinatorial_map<3, 3,...>
    The dimension of the combinatorial map and the dimension of the ambient space have to be 3.

    \param lcc
    an instance of LCC

    \pre `successful partitioning`
   */
  template<typename LCC>
  void get_linear_cell_complex(LCC& lcc) const;

  /// @}
}
#else

template<typename GeomTraits, typename Intersection_Kernel>
class Kinetic_shape_reconstruction_3 {

public:
  using Kernel = GeomTraits;
  using LCC = Linear_cell_complex_for_combinatorial_map<3, 3>;

private:
  using FT      = typename Kernel::FT;
  using Point_3 = typename Kernel::Point_3;

  using Data_structure = KSR_3::Data_structure<Kernel, Intersection_Kernel>;

  using IVertex = typename Data_structure::IVertex;
  using IEdge   = typename Data_structure::IEdge;

  using From_exact = CGAL::Cartesian_converter<Intersection_Kernel, Kernel>;

  using Initializer = KSR_3::Initializer<Kernel, Intersection_Kernel>;
  using Propagation = KSR_3::FacePropagation<Kernel, Intersection_Kernel>;
  using Finalizer   = KSR_3::Finalizer<Kernel, Intersection_Kernel>;

  using Polygon_mesh = CGAL::Surface_mesh<Point_3>;
  using Vertex_index = typename Polygon_mesh::Vertex_index;
  using Timer        = CGAL::Real_timer;
  using Parameters   = KSR::Parameters_3<FT>;

private:
  Parameters m_parameters;
  Data_structure m_data;
  std::size_t m_num_events;

public:
  Kinetic_shape_reconstruction_3(
    const bool verbose = true,
    const bool debug   = false) :
  m_parameters(verbose, debug, false), // use true here to export all steps
  m_data(m_parameters),
  m_num_events(0)
  { }

  template<
    typename InputRange,
    typename PolygonMap,
    typename NamedParameters>
  bool partition(
    const InputRange & input_range,
    const PolygonMap polygon_map,
    const NamedParameters & np) {

    Timer timer;
    m_parameters.k = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::k_intersections), 1);
    m_parameters.n = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::n_subdivisions), 0);
    m_parameters.enlarge_bbox_ratio = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::enlarge_bbox_ratio), FT(11) / FT(10));
    m_parameters.distance_tolerance = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::distance_tolerance), FT(5) / FT(10));
    m_parameters.reorient = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::reorient), false);

    std::cout.precision(20);
    if (input_range.size() == 0) {
      CGAL_warning_msg(input_range.size() > 0,
        "WARNING: YOUR INPUT IS EMPTY! RETURN WITH NO CHANGE!");
      return false;
    }

    if (m_parameters.n != 0) {
      CGAL_assertion_msg(false, "TODO: IMPLEMENT KINETIC SUBDIVISION!");
      if (m_parameters.n > 3) {
        CGAL_warning_msg(m_parameters.n <= 3,
          "WARNING: DOES IT MAKE SENSE TO HAVE MORE THAN 64 INPUT BLOCKS? SETTING N TO 3!");
        m_parameters.n = 3;
      }
    }

    if (m_parameters.enlarge_bbox_ratio < FT(1)) {
      CGAL_warning_msg(m_parameters.enlarge_bbox_ratio >= FT(1),
        "WARNING: YOU SET ENLARGE_BBOX_RATIO < 1.0! THE VALID RANGE IS [1.0, +INF). SETTING TO 1.0!");
      m_parameters.enlarge_bbox_ratio = FT(1);
    }

    if (m_parameters.verbose) {
      const unsigned int num_blocks = static_cast<unsigned int>(std::pow(m_parameters.n + 1, 3));
      const std::string is_reorient = (m_parameters.reorient ? "true" : "false");

      std::cout << std::endl << "--- PARTITION OPTIONS: " << std::endl;
      std::cout << "* number of intersections k: " << m_parameters.k << std::endl;
      std::cout << "* number of subdivisions per bbox side: " << m_parameters.n << std::endl;
      std::cout << "* number of subdivision blocks: " << num_blocks << std::endl;
      std::cout << "* enlarge bbox ratio: " << m_parameters.enlarge_bbox_ratio << std::endl;
      std::cout << "* reorient: " << is_reorient << std::endl;
    }

    if (m_parameters.verbose) {
      std::cout << std::endl << "--- INITIALIZING PARTITION:" << std::endl;
    }

    // Initialization.
    timer.reset();
    timer.start();
    m_data.clear();

    Initializer initializer(m_data, m_parameters);
    initializer.initialize(input_range, polygon_map);

    timer.stop();
    const double time_to_initialize = timer.time();

    // if (m_parameters.verbose) {
    //   std::cout << std::endl << "* initialization (sec.): " << time_to_initialize << std::endl;
    //   std::cout << "INITIALIZATION SUCCESS!" << std::endl << std::endl;
    // }
    // exit(EXIT_SUCCESS);

    // Output planes.
    // for (std::size_t i = 6; i < m_data.number_of_support_planes(); ++i) {
    //   std::cout << m_data.support_plane(i).plane() << std::endl;
    // }

    if (m_parameters.k == 0) { // for k = 0, we skip propagation
      CGAL_warning_msg(m_parameters.k > 0,
        "WARNING: YOU SET K TO 0! THAT MEANS NO PROPAGATION! THE VALID VALUES ARE {1,2,...}. INTERSECT AND RETURN!");
      return false;
    }

    if (m_parameters.verbose) {
      std::cout << std::endl << "--- RUNNING THE QUEUE:" << std::endl;
      std::cout << "* propagation started" << std::endl;
    }

    // Propagation.
    timer.reset();
    timer.start();
    std::size_t num_queue_calls = 0;

    Propagation propagation(m_data, m_parameters);
    std::tie(num_queue_calls, m_num_events) = propagation.propagate();

    timer.stop();
    const double time_to_propagate = timer.time();

    if (m_parameters.verbose) {
      std::cout << "* propagation finished" << std::endl;
      std::cout << "* number of queue calls: " << num_queue_calls << std::endl;
      std::cout << "* number of events handled: " << m_num_events << std::endl;
    }

    if (m_parameters.verbose) {
      std::cout << std::endl << "--- FINALIZING PARTITION:" << std::endl;
    }

    // Finalization.
    timer.reset();
    timer.start();
    if (m_parameters.debug)
      dump(m_data, "final-" + m_parameters.k);

    Finalizer finalizer(m_data, m_parameters);

    if (m_parameters.verbose)
      std::cout << "* checking final mesh integrity ...";

    CGAL_assertion(m_data.check_integrity(true, true, true));

    if (m_parameters.verbose)
      std::cout << " done" << std::endl;

    if (m_parameters.verbose)
      std::cout << "* getting volumes ..." << std::endl;

    finalizer.create_polyhedra();
    timer.stop();
    const double time_to_finalize = timer.time();

    if (m_parameters.verbose) {
      std::cout << "* found all together " << m_data.number_of_volumes() << " volumes" << std::endl;

      for (std::size_t i = 0; i < m_data.number_of_support_planes(); i++) {
        dump_2d_surface_mesh(m_data, i, "final-surface-mesh-" + std::to_string(i));
      }
    }

    // Timing.
    if (m_parameters.verbose) {
      std::cout << std::endl << "--- TIMING (sec.):" << std::endl;
    }
    const double total_time =
      time_to_initialize + time_to_propagate + time_to_finalize;

    std::cout << "* initialization: " << time_to_initialize << std::endl;
    std::cout << "* propagation: " << time_to_propagate << std::endl;
    std::cout << "* finalization: " << time_to_finalize << std::endl;
    std::cout << "* total time: " << total_time << std::endl;

    return true;
  }

  template<
    typename InputRange,
    typename PointMap,
    typename VectorMap,
    typename SemanticMap,
    typename NamedParameters>
  bool reconstruct(
    const InputRange& input_range,
    const PointMap point_map,
    const VectorMap normal_map,
    const SemanticMap semantic_map,
    const std::string file_name,
    const NamedParameters& np) {

    using Reconstruction = KSR_3::Reconstruction<
      InputRange, PointMap, VectorMap, SemanticMap, Kernel, Intersection_Kernel>;

    Reconstruction reconstruction(
      input_range, point_map, normal_map, semantic_map, m_data, m_parameters.verbose, m_parameters.debug);

    bool success = reconstruction.detect_planar_shapes(file_name, np);
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

  template<
    typename InputRange,
    typename PointMap,
    typename VectorMap,
    typename SemanticMap,
    typename RegionMap,
    typename NamedParameters>
  bool reconstruct(
    const InputRange& input_range,
    const PointMap point_map,
    const VectorMap normal_map,
    const SemanticMap semantic_map,
    const RegionMap region_map,
    const NamedParameters& np) {

    using Reconstruction = KSR_3::Reconstruction<
      InputRange, PointMap, VectorMap, SemanticMap, Kernel, Intersection_Kernel>;

    Reconstruction reconstruction(
      input_range, point_map, normal_map, semantic_map, m_data, m_parameters.verbose, m_parameters.debug);

    bool success = reconstruction.planar_shapes_from_map(region_map, np);
    if (!success) {
      CGAL_assertion_msg(false, "ERROR: RECONSTRUCTION, DETECTING PLANAR SHAPES FAILED!");
      return false;
    }
    // exit(EXIT_SUCCESS);

    //success = reconstruction.regularize_planar_shapes(np);
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

  std::size_t number_of_events() const {
    return m_num_events;
  }

  int number_of_support_planes() const {
    return static_cast<int>(m_data.number_of_support_planes());
  }

  std::size_t number_of_vertices(const int support_plane_idx = -1) const {

    CGAL_assertion(support_plane_idx < number_of_support_planes());
    if (support_plane_idx >= number_of_support_planes()) return std::size_t(-1);
    if (support_plane_idx < 0) {
      return m_data.igraph().number_of_vertices();
    }

    CGAL_assertion(support_plane_idx >= 0);
    const std::size_t sp_idx = static_cast<std::size_t>(support_plane_idx);
    return static_cast<std::size_t>(m_data.mesh(sp_idx).number_of_vertices());
  }

  std::size_t number_of_edges(const int support_plane_idx = -1) const {

    CGAL_assertion(support_plane_idx < number_of_support_planes());
    if (support_plane_idx >= number_of_support_planes()) return std::size_t(-1);
    if (support_plane_idx < 0) {
      return m_data.igraph().number_of_edges();
    }

    CGAL_assertion(support_plane_idx >= 0);
    const std::size_t sp_idx = static_cast<std::size_t>(support_plane_idx);
    return static_cast<std::size_t>(m_data.mesh(sp_idx).number_of_edges());
  }

  std::size_t number_of_faces(const int support_plane_idx = -1) const {

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

  std::size_t number_of_volumes() const {
    return m_data.volumes().size();
  }

  int support_plane_index(const std::size_t polygon_index) const {
    const int support_plane_idx = m_data.support_plane_index(polygon_index);
    CGAL_assertion(support_plane_idx >= 6);
    return support_plane_idx;
  }

  template<typename LCC>
  void get_linear_cell_complex(LCC &lcc) const {
    lcc.clear();

    From_exact from_exact;
    std::vector<bool> used_vertices(m_data.igraph().number_of_vertices(), false);
    std::vector<int> remap(m_data.igraph().number_of_vertices(), -1);
    std::vector<Point_3> mapped_vertices;
    mapped_vertices.reserve(m_data.igraph().number_of_vertices());

    for (const auto& volume : m_data.volumes()) {
      for (const auto& vertex : volume.pvertices) {
        CGAL_assertion(m_data.has_ivertex(vertex));
        IVertex ivertex = m_data.ivertex(vertex);
        if (remap[ivertex] == -1) {
          remap[ivertex] = static_cast<int>(mapped_vertices.size());
          mapped_vertices.push_back(from_exact(m_data.point_3(ivertex)));
        }
      }
    }

    CGAL::Linear_cell_complex_incremental_builder_3<LCC> ib(lcc);
    for (const auto& p : mapped_vertices)
      ib.add_vertex(p);

    for (const auto& vol : m_data.volumes()) {
      ib.begin_surface();
      for (std::size_t i = 0; i < vol.pfaces.size(); i++) {
        auto vertex_range = m_data.pvertices_of_pface(vol.pfaces[i]);
        ib.begin_facet();
        if (vol.pface_oriented_outwards[i]) {
          typename Data_structure::PVertex_of_pface_iterator it = vertex_range.begin();
          while (it != vertex_range.end()) {
            CGAL_assertion(m_data.has_ivertex(*it));
            IVertex ivertex = m_data.ivertex(*it);
            ib.add_vertex_to_facet(static_cast<std::size_t>(remap[ivertex]));
            it++;
          }
        }
        else {
          typename Data_structure::PVertex_of_pface_iterator it = vertex_range.end()--;
          do {
            CGAL_assertion(m_data.has_ivertex(*it));
            IVertex ivertex = m_data.ivertex(*it);
            ib.add_vertex_to_facet(static_cast<std::size_t>(remap[ivertex]));
            if (it == vertex_range.begin())
              break;
            it--;
          } while (true);
        }
        ib.end_facet();
      }
      ib.end_surface();
    }

    lcc.display_characteristics(std::cout) << std::endl;
  }

  /*******************************
  **           OUTPUT           **
  ********************************/

  template<typename VertexOutputIterator>
  VertexOutputIterator output_partition_vertices(
    VertexOutputIterator vertices, const int support_plane_idx = -1) const {
    From_exact from_EK;

    CGAL_assertion(support_plane_idx < number_of_support_planes());
    if (support_plane_idx >= number_of_support_planes()) return vertices;
    if (support_plane_idx < 0) {
      const auto all_ivertices = m_data.ivertices();
      for (const auto ivertex : all_ivertices) {
        *(vertices++) = from_EK(m_data.point_3(ivertex));
      }
      return vertices;
    }

    CGAL_assertion(support_plane_idx >= 0);
    const std::size_t sp_idx = static_cast<std::size_t>(support_plane_idx);
    const auto all_pvertices = m_data.pvertices(sp_idx);
    for (const auto pvertex : all_pvertices) {
      CGAL_assertion(m_data.has_ivertex(pvertex));
      const auto ivertex = m_data.ivertex(pvertex);
      *(vertices++) = from_EK(m_data.point_3(ivertex));
    }
    return vertices;
  }

  template<typename EdgeOutputIterator>
  EdgeOutputIterator output_partition_edges(
    EdgeOutputIterator edges, const int support_plane_idx = -1) const {
    From_exact from_EK;

    CGAL_assertion(support_plane_idx < number_of_support_planes());
    if (support_plane_idx >= number_of_support_planes()) return edges;
    if (support_plane_idx < 0) {
      const auto all_iedges = m_data.iedges();
      for (const auto iedge : all_iedges) {
        *(edges++) = from_EK(m_data.segment_3(iedge));
      }
      return edges;
    }

    CGAL_assertion(support_plane_idx >= 0);
    const std::size_t sp_idx = static_cast<std::size_t>(support_plane_idx);
    const auto all_pedges = m_data.pedges(sp_idx);
    for (const auto pedge : all_pedges) {
      CGAL_assertion(m_data.has_iedge(pedge));
      const auto iedge = m_data.iedge(pedge);
      *(edges++) = from_EK(m_data.segment_3(iedge));
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
    From_exact from_EK;

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
        polygon_mesh.add_vertex(from_EK(m_data.point_3(ivertex)));
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
    VolumeOutputIterator volumes) const {
    for (std::size_t i = 0; i < m_data.number_of_volumes(); ++i) {
      output_partition_volume(volumes, i);
    }
    return volumes;
  }

  template<typename VolumeOutputIterator>
  VolumeOutputIterator output_partition_volume(
    VolumeOutputIterator volumes, const std::size_t volume_index) const {

    CGAL_assertion(volume_index < number_of_volumes());
    if (volume_index >= number_of_volumes()) return volumes;

    std::vector<EK::Point_3> vertices;
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

    CGAL_assertion(volume_index < number_of_volumes());
    if (volume_index >= number_of_volumes()) return;

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

  template<typename VertexOutputIterator, typename FaceOutputIterator>
  void output_reconstructed_model(
    VertexOutputIterator vertices, FaceOutputIterator faces) const {
    From_exact from_EK;

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
          *(vertices++) = from_EK(m_data.point_3(ivertex));
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
};

#endif //DOXYGEN_RUNNING

} // namespace CGAL

#endif // CGAL_KINETIC_SHAPE_RECONSTRUCTION_3_H
