// Copyright (c) 2023 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sven Oesau, Florent Lafarge, Dmitry Anisimov, Simon Giraudot

#ifndef CGAL_KINETIC_SHAPE_PARTITION_3_H
#define CGAL_KINETIC_SHAPE_PARTITION_3_H

// #include <CGAL/license/Kinetic_shape_reconstruction.h>

// Boost includes.
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Named_function_parameters.h>

// CGAL includes.
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
//#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
//#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
//#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Real_timer.h>

#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_incremental_builder.h>
//#include <CGAL/draw_linear_cell_complex.h>

// Internal includes.
#include <CGAL/KSR/utils.h>
#include <CGAL/KSR/debug.h>
#include <CGAL/KSR/parameters.h>

#include <CGAL/KSR_3/Data_structure.h>
//#include <CGAL/KSR_3/Reconstruction.h>
#include <CGAL/KSR_3/Initializer.h>
#include <CGAL/KSR_3/FacePropagation.h>
#include <CGAL/KSR_3/Finalizer.h>

namespace CGAL {

/*!
* \ingroup PkgKineticShapeReconstructionRef
  \brief creates the kinetic partition of the bounding box of the polygons given as input data. Use `Kinetic_shape_partition_3()` to create an empty object, `insert()` to provide input data and `initialize()` to prepare the partition or use `Kinetic_shape_partition_3()`

  \tparam GeomTraits
    must be a model of `KineticShapePartitionTraits_3`.

  \tparam IntersectionTraits
    must be a model of `Kernel` using exact computations. Defaults to `CGAL::Exact_predicates_exact_constructions_kernel`.
*/
template<typename GeomTraits, typename IntersectionTraits = CGAL::Exact_predicates_exact_constructions_kernel>
class Kinetic_shape_partition_3 {

public:
  using Kernel = typename GeomTraits;
  using Intersection_kernel = IntersectionTraits;

  //using Point_3 = typename GeomTraits::Point_3;

private:
  using FT      = typename GeomTraits::FT;

  using Data_structure = KSR_3::Data_structure<Traits>;

  using IVertex = typename Data_structure::IVertex;
  using IEdge   = typename Data_structure::IEdge;

  using From_exact = typename CGAL::Cartesian_converter<Intersection_kernel, Kernel>;

  using Initializer = KSR_3::Initializer<Kernel, Intersection_kernel>;
  using Propagation = KSR_3::FacePropagation<Kernel, Intersection_kernel>;
  using Finalizer   = KSR_3::Finalizer<Kernel, Intersection_kernel>;

  using Polygon_mesh = CGAL::Surface_mesh<Point_3>;
  using Timer        = CGAL::Real_timer;
  using Parameters   = KSR::Parameters_3<FT>;

private:
  Parameters m_parameters;
  Data_structure m_data;
  std::size_t m_num_events;

public:
  /// \name Initialization
  /// @{
  /*!
  \brief constructs an empty kinetic shape partition object. Use `insert()` afterwards to insert polygons into the partition and `initialize()` to initialize the partition.

  \tparam NamedParameters
  a sequence of \ref bgl_namedparameters "Named Parameters"

  \param np
  a sequence of \ref bgl_namedparameters "Named Parameters"
  among the ones listed below

  \cgalNamedParamsBegin
    \cgalParamNBegin{verbose}
      \cgalParamDescription{Write timing and internal information to std::cout.}
      \cgalParamType{bool}
      \cgalParamDefault{false}
    \cgalParamNEnd
    \cgalParamNBegin{debug}
      \cgalParamDescription{Export of intermediate results.}
      \cgalParamType{bool}
      \cgalParamDefault{false}
    \cgalParamNEnd
  \cgalNamedParamsEnd
  */
  template<typename NamedParameters = parameters::Default_named_parameters>
  Kinetic_shape_partition_3(
    const NamedParameters& np = CGAL::parameters::default_values()) :
    m_parameters(np, false), // use true here to export all steps
    m_data(m_parameters),
    m_num_events(0)
  { }

  /*!
  \brief constructs a kinetic shape partition object and initializes it.

  \tparam InputRange
  must be a model of `ConstRange` whose iterator type is `RandomAccessIterator` and whose value type is Point_3.

  \tparam PolygonRange
  contains index ranges to form polygons by providing indices into InputRange

  \tparam NamedParameters
  a sequence of \ref bgl_namedparameters "Named Parameters"

  \param input_range
  an instance of `InputRange` with 3D points and corresponding 3D normal vectors

  \param polygon_range
  a range of polygons defined by a range of indices into `input_range`

  \param np
  a sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

  \pre !input_range.empty() and !polygon_map.empty()

  \cgalNamedParamsBegin
     \cgalParamNBegin{point_map}
       \cgalParamDescription{a property map associating points to the elements of the `input_range`}
       \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type of the iterator of `InputRange` and whose value type is `GeomTraits::Point_3`}
       \cgalParamDefault{`CGAL::Identity_property_map<GeomTraits::Point_3>`}
     \cgalParamNEnd
    \cgalParamNBegin{verbose}
      \cgalParamDescription{Write timing and internal information to std::cout.}
      \cgalParamType{bool}
      \cgalParamDefault{false}
    \cgalParamNEnd
    \cgalParamNBegin{reorient_bbox}
      \cgalParamDescription{Use the oriented bounding box instead of the axis-aligned bounding box. While the z direction is maintained, the x axis is aligned with the largest variation in the horizontal plane.}
      \cgalParamType{bool}
      \cgalParamDefault{false}
    \cgalParamNEnd
    \cgalParamNBegin{bbox_dilation_ratio}
      \cgalParamDescription{Factor for extension of the bounding box of the input data to be used for the partition.}
      \cgalParamType{FT}
      \cgalParamDefault{1.1}
    \cgalParamNEnd
    \cgalParamNBegin{angle_tolerance}
      \cgalParamDescription{The tolerance angle to snap the planes of two input polygons into one plane.}
      \cgalParamType{FT}
      \cgalParamDefault{5}
    \cgalParamNEnd
    \cgalParamNBegin{distance_tolerance}
      \cgalParamDescription{The tolerance distance to snap the planes of two input polygons into one plane.}
      \cgalParamType{FT}
      \cgalParamDefault{1% of the bounding box diagonal}
    \cgalParamNEnd
  \cgalNamedParamsEnd
  */
  template<
    typename InputRange,
    typename PolygonRange,
    typename NamedParameters = parameters::Default_named_parameters>
  Kinetic_shape_partition_3(
    const InputRange& input_range,
    const PolygonRange polygon_range,
    const NamedParameters & np = CGAL::parameters::default_values()) :
    m_parameters(np, false), // use true here to export all steps
    m_data(m_parameters),
    m_num_events(0)
  {
    initialize(input_range, polygon_range, np);
  }

  /*!
  \brief inserts polygons, requires initialize() afterwards to have effect.

  \tparam InputRange
  must be a model of `ConstRange` whose iterator type is `RandomAccessIterator` and whose value type is Point_3.

  \tparam PolygonRange
  contains index ranges to form polygons by providing indices into InputRange

  \tparam NamedParameters
  a sequence of \ref bgl_namedparameters "Named Parameters"

  \param input_range
  an instance of `InputRange` with 3D points

  \param polygon_range
  a range of polygons defined by a range of indices into `input_range`

  \param np
  a sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

  \cgalNamedParamsBegin
     \cgalParamNBegin{point_map}
       \cgalParamDescription{a property map associating points to the elements of the `input_range`}
       \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type of the iterator of `InputRange` and whose value type is `GeomTraits::Point_3`}
       \cgalParamDefault{`CGAL::Identity_property_map<GeomTraits::Point_3>`}
     \cgalParamNEnd
  \cgalNamedParamsEnd
  */

  template<typename InputRange, typename PolygonRange, typename NamedParameters = parameters::Default_named_parameters>
  void insert(
    const InputRange& input_range,
    const PolygonRange polygon_range,
    const NamedParameters& np = CGAL::parameters::default_values()) {}

  /*!
  \brief initializes the kinetic partition of the bounding box.

  \tparam NamedParameters
  a sequence of \ref bgl_namedparameters "Named Parameters"

  \param np
  a sequence of \ref bgl_namedparameters "Named Parameters"
  among the ones listed below

  \pre Input data has been provided via `insert()`.

  \cgalNamedParamsBegin
    \cgalParamNBegin{reorient_bbox}
      \cgalParamDescription{Use the oriented bounding box instead of the axis-aligned bounding box.}
      \cgalParamType{bool}
      \cgalParamDefault{false}
    \cgalParamNEnd
    \cgalParamNBegin{bbox_dilation_ratio}
      \cgalParamDescription{Factor for extension of the bounding box of the input data to be used for the partition.}
      \cgalParamType{FT}
      \cgalParamDefault{1.1}
    \cgalParamNEnd
    \cgalParamNBegin{angle_tolerance}
      \cgalParamDescription{The tolerance angle to snap the planes of two input polygons into one plane.}
      \cgalParamType{FT}
      \cgalParamDefault{5}
    \cgalParamNEnd
    \cgalParamNBegin{distance_tolerance}
      \cgalParamDescription{The tolerance distance to snap the planes of two input polygons into one plane.}
      \cgalParamType{FT}
      \cgalParamDefault{0.5}
    \cgalParamNEnd
  \cgalNamedParamsEnd
  */

  template<
    typename NamedParameters = parameters::Default_named_parameters>
  void initialize(
    const NamedParameters& np = CGAL::parameters::default_values()) {

    Timer timer;
    m_parameters.bbox_dilation_ratio = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::bbox_dilation_ratio), FT(11) / FT(10));
    m_parameters.angle_tolerance = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::angle_tolerance), FT(5) / FT(10));
    m_parameters.distance_tolerance = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::distance_tolerance), FT(5) / FT(10));
    m_parameters.reorient_bbox = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::reorient_bbox), false);

    std::cout.precision(20);
    if (input_range.size() == 0) {
      CGAL_warning_msg(input_range.size() > 0,
        "Warning: Your input is empty!");
      return false;
    }

    if (m_parameters.bbox_dilation_ratio < FT(1)) {
      CGAL_warning_msg(m_parameters.bbox_dilation_ratio >= FT(1),
        "Warning: You set enlarge_bbox_ratio < 1.0! The valid range is [1.0, +inf). Setting to 1.0!");
      m_parameters.bbox_dilation_ratio = FT(1);
    }

    if (m_parameters.verbose) {
      const unsigned int num_blocks = static_cast<unsigned int>(std::pow(m_parameters.n + 1, 3));
      //const std::string is_reorient = (m_parameters.reorient ? "true" : "false");

      std::cout << std::endl << "--- PARTITION OPTIONS: " << std::endl;
      std::cout << "* enlarge bbox ratio: " << m_parameters.bbox_dilation_ratio << std::endl;
    }

    if (m_parameters.verbose) {
      std::cout << std::endl << "--- INITIALIZING PARTITION:" << std::endl;

      // Initialization.
      timer.reset();
      timer.start();
    }

    m_data.clear();

    Initializer initializer(m_data, m_parameters);
    initializer.initialize(input_range, polygon_map);

    // Timing.
    if (m_parameters.verbose) {
      timer.stop();
      const double time_to_initialize = timer.time();
      std::cout << "* initialization time: " << time_to_initialize << std::endl;
    }

    return true;
  }

  /*!
  \brief propagates the kinetic polygons in the initialized partition.

  \param k
   maximum number of allowed intersections for each input polygon before its expansion stops.

  \pre successful initialization and k != 0
  */
  void partition(std::size_t k) {
    Timer timer;
    std::cout.precision(20);

    // Already initialized?
    if (m_data.number_of_support_planes() < 6) {
      std::cout << "Kinetic partition not initialized or empty. Number of support planes: " << m_data.number_of_support_planes() << std::endl;

      return;
    }

    if (k == 0) { // for k = 0, we skip propagation
      std::cout << "k needs to be a positive number" << std::endl;

      return;
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
    std::tie(num_queue_calls, m_num_events) = propagation.propagate(k);

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
      dump(m_data, "final-" + std::to_string(m_parameters.k));

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

    if (m_parameters.verbose)
      std::cout << "* found all together " << m_data.number_of_volumes() << " volumes" << std::endl;

    if (m_parameters.debug)
      for (std::size_t i = 0; i < m_data.number_of_support_planes(); i++)
        dump_2d_surface_mesh(m_data, i, "final-surface-mesh-" + std::to_string(i));

    // Timing.
    if (m_parameters.verbose) {
      std::cout << std::endl << "--- TIMING (sec.):" << std::endl;

      std::cout << "* propagation: " << time_to_propagate << std::endl;
      std::cout << "* finalization: " << time_to_finalize << std::endl;
    }

    return;
  }

  /// @}

  /*******************************
  **         Access         **
  ********************************/

  /// \name Access
  /// @{
  /*!
  \brief returns the number of support planes in the kinetic partition. They originate from the planes of the input polygons and the bounding box.

  \pre successful partition
  */

  int number_of_support_planes() const {
    return static_cast<int>(m_data.number_of_support_planes());
  }

  /*!
  \brief returns the number of vertices in the kinetic partition.

  \pre successful partition
  */
  std::size_t number_of_vertices() const {
    return 0;
  }

  /*!
  \brief returns the number of faces in the kinetic partition.

  \pre successful partition
  */
  std::size_t number_of_faces() const {
    return 0;
  }

  /*!
  \brief returns the number of volumes created by the kinetic partition.

  \pre successful partition
  */
  std::size_t number_of_volumes() const {
    return m_data.volumes().size();
  }

#ifndef DOXYGEN_RUNNING
  /*!
  \brief Point vector for mapping vertex indices to positions.

  @return
   vector of points.

    \pre successful partition
  */
  const std::vector<Point_3>& vertices() const;

  /*!
  \brief Vertex indices of face.

  \param face_index
   index of the query face.

  @return
   vector of vertex indices.

  \pre successful partition
  */
  const std::vector<std::size_t>& vertices(std::size_t face_index) const;

  /*!
  \brief Face indices of the volume.

  \param volume_index
   index of the query volume.

  @return
   vector of face indices.

  \pre successful partition
  */
  const std::vector<std::size_t>& face(std::size_t volume_index) const {
    CGAL_assertion(m_data.number_of_volumes() > volume_index);
    return m_data.volumes()[volume_index].faces;
  }

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

    \pre successful partition
  */
  const std::pair<int, int>& neighbors(std::size_t face_index) const {
    CGAL_assertion(m_data.number_of_volumes() > volume_index);
    return m_data.
  }

  /*!
  \brief Retrieves the support plane generated from the input polygon.

  \param input_polygon_index
   index of the input polygon.

  @return
   index into polygon_map provided on initialization.

  \pre successful partition
  */
  std::size_t support_plane_index(const std::size_t input_polygon_index) const {
      const int support_plane_idx = m_data.support_plane_index(input_polygon_index);
      CGAL_assertion(support_plane_idx >= 6);
      return support_plane_idx;
  }

#endif

  /*!
   \brief creates a linear cell complex from the kinetic partition.

   \tparam LCC
    most be a model of `LinearCellComplex`
    The dimension of the combinatorial map and the dimension of the ambient space have to be 3.

   \param lcc
    an instance of LCC

   \pre successful partition
  */

  template<typename LCC>
  void get_linear_cell_complex(LCC& lcc) const {
    using LCC_Kernel = typename LCC::Traits;
    CGAL::Cartesian_converter<Intersection_kernel, LCC_Kernel> conv;
    lcc.clear();

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
          mapped_vertices.push_back(conv(m_data.point_3(ivertex)));
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
            it--;
            if (it == vertex_range.begin())
              break;
          } while (true);
        }
        ib.end_facet();
      }
      ib.end_surface();
    }

    lcc.display_characteristics(std::cout) << std::endl;
  }

  /// @}

  /*
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
    From_exact from_EK;

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
          *(vertices++) = from_EK(m_data.point_3(ivertex));
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
*/

  /*******************************
  **           MEMORY           **
  ********************************/
  /*!
   \brief clears all input data and the kinetic partition.
   */
  void clear() {
    m_data.clear();
    m_num_events = 0;
  }

private:
  /*

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
  }*/
};

} // namespace CGAL

#endif // CGAL_KINETIC_SHAPE_PARTITION_3_H
