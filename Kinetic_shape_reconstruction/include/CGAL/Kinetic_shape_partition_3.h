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

#include <numeric>

// CGAL includes.
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Real_timer.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Base_with_time_stamp.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>

#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_incremental_builder_3.h>
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

//#define OVERLAY_2_DEBUG
#define OVERLAY_2_CHECK

#include "cdtLower.h"

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
  using Kernel = GeomTraits;
  using Intersection_kernel = IntersectionTraits;

  using Point_3 = typename Kernel::Point_3;

private:
  using FT      = typename Kernel::FT;
  using Point_2 = typename Kernel::Point_2;
  using Plane_3 = typename Kernel::Plane_3;
  using Line_3 = typename Kernel::Line_3;
  using Line_2 = typename Kernel::Line_2;
  using Transform_3 = CGAL::Aff_transformation_3<Kernel>;

  using Data_structure = KSR_3::Data_structure<Kernel, Intersection_kernel>;

  using IVertex = typename Data_structure::IVertex;
  using IEdge   = typename Data_structure::IEdge;

  using From_exact = typename CGAL::Cartesian_converter<Intersection_kernel, Kernel>;
  using To_exact = typename CGAL::Cartesian_converter<Kernel, Intersection_kernel>;

  using Initializer = KSR_3::Initializer<Kernel, Intersection_kernel>;
  using Propagation = KSR_3::FacePropagation<Kernel, Intersection_kernel>;
  using Finalizer   = KSR_3::Finalizer<Kernel, Intersection_kernel>;

  using Polygon_mesh = CGAL::Surface_mesh<Point_3>;
  using Timer        = CGAL::Real_timer;
  using Parameters   = KSR::Parameters_3<FT>;

  struct VI
  {
    VI()
      : input(false), idx(-1)
    {}

    void set_index(std::size_t i) {
      idx = i;
    }

    void set_point(const typename Intersection_kernel::Point_3& p)
    {
      point_3 = p;
      input = true;
    }

    typename Intersection_kernel::Point_3 point_3;
    std::size_t idx;
    bool input;
  };

  // Each face gets as id
  // The overlay face gets also the id from A and B
  struct ID {
    ID()
      : id(-1), idA(-1), idB(-1)
    {}

    ID(const ID& other)
      :id(other.id), idA(other.idA), idB(other.idB)
    {}

    ID& operator=(const ID& other)
    {
      id = other.id;
      idA = other.idA;
      idB = other.idB;
      return *this;
    }

    int id, idA, idB;
  };

  typedef CGAL::Triangulation_vertex_base_with_info_2<VI, Intersection_kernel> Vbi;
  typedef CGAL::Triangulation_face_base_with_info_2<ID, Intersection_kernel> Fbi;
  typedef CGAL::Constrained_triangulation_face_base_2<Intersection_kernel, Fbi>    Fb;
  typedef CGAL::Triangulation_data_structure_2<Vbi, Fb>       TDS;
  typedef CGAL::Exact_intersections_tag                     Itag;
  typedef CGAL::Constrained_Delaunay_triangulation_2<Intersection_kernel, TDS, Itag> CDT;
  typedef CGAL::Constrained_triangulation_plus_2<CDT>       CDTplus;
  typedef typename CDTplus::Vertex_handle                            Vertex_handle;
  typedef typename CDTplus::Finite_vertices_iterator                 Finite_vertices_iterator;
  typedef typename CDTplus::Finite_faces_iterator                    Finite_faces_iterator;

private:
  struct Sub_partition {
    Sub_partition() : parent(-1) {}
    std::shared_ptr<Data_structure> m_data;
    std::array<Point_3, 8> bbox;
    std::vector<std::size_t> input_polygons;
    std::size_t parent;
    std::vector<std::size_t> children;
    std::size_t split_plane;

    std::vector<typename Data_structure::Volume_cell> volumes;
    std::vector<std::vector<std::size_t> > face2vertices;
    std::vector<Point_3> exact_vertices;


    // Merged data from children
  };

  Parameters m_parameters;
  std::array<Point_3, 8> m_bbox;
  std::vector<Sub_partition> m_partition_nodes; // Tree of partitions.
  std::vector<std::size_t> m_partitions; // Contains the indices of the leaf nodes, the actual partitions to be calculated.
  std::size_t m_num_events;
  std::vector<std::vector<Point_3> > m_input_polygons;

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
    m_parameters(
      parameters::choose_parameter(parameters::get_parameter(np, internal_np::verbose), false),
      parameters::choose_parameter(parameters::get_parameter(np, internal_np::debug), false)), // use true here to export all steps
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
    m_parameters(
      parameters::choose_parameter(parameters::get_parameter(np, internal_np::verbose), false),
      parameters::choose_parameter(parameters::get_parameter(np, internal_np::debug), false)), // use true here to export all steps
    m_data(m_parameters),
    m_num_events(0)
  {
    insert(input_range, polygon_range, np);
    initialize(np);
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
    const NamedParameters& np = CGAL::parameters::default_values()) {
    for (auto poly : polygon_range) {
      std::vector<Point_3> pts;
      pts.reserve(poly.size());
      for (auto it : poly)
        pts.push_back(*(input_range.begin() + it));

      Plane_3 pl;
      CGAL::linear_least_squares_fitting_3(pts.begin(), pts.end(), pl, CGAL::Dimension_tag<0>());

      std::vector<Point_2> pts2d(pts.size());
      for (std::size_t i = 0; i < pts.size(); i++)
        pts2d[i] = pl.to_2d(pts[i]);

      std::vector<Point_2> ch;
      CGAL::convex_hull_2(pts2d.begin(), pts2d.end(), std::back_inserter(ch));

      m_input_polygons.push_back(std::vector<Point_3>(ch.size()));

      for (std::size_t i = 0; i < ch.size(); i++)
        m_input_polygons.back()[i] = pl.to_3d(ch[i]);
    }
  }

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
    if (m_input_polygons.size() == 0) {
      std::cout << "Warning: Your input is empty!";
      return;
    }

    if (m_parameters.bbox_dilation_ratio < FT(1)) {
      CGAL_warning_msg(m_parameters.bbox_dilation_ratio >= FT(1),
        "Warning: You set enlarge_bbox_ratio < 1.0! The valid range is [1.0, +inf). Setting to 1.0!");
      m_parameters.bbox_dilation_ratio = FT(1);
    }

    if (m_parameters.verbose) {
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

    m_partition_nodes.resize(1);
    create_bounding_box(m_parameters.bbox_dilation_ratio, m_parameters.reorient_bbox, m_partition_nodes[0].bbox);

    m_partition_nodes[0].input_polygons.resize(m_input_polygons.size());
    std::iota(m_partition_nodes[0].input_polygons.begin(), m_partition_nodes[0].input_polygons.end(), 0);

    split_partition(0);

    for (std::size_t idx : m_partitions) {
      Sub_partition& partition = m_partition_nodes[idx];

      partition.m_data = std::make_shared<Data_structure>(m_parameters, std::to_string(idx) + "-");

      std::vector<std::vector<Point_3> > input_polygons(partition.input_polygons.size());
      for (std::size_t i = 0; i < partition.input_polygons.size(); i++)
        input_polygons[i] = m_input_polygons[partition.input_polygons[i]];

      Initializer initializer(input_polygons, *partition.m_data, m_parameters);
      initializer.initialize(partition.bbox);
    }

    // Timing.
    if (m_parameters.verbose) {
      timer.stop();
      const double time_to_initialize = timer.time();
      std::cout << "* initialization time: " << time_to_initialize << std::endl;
    }
  }

  /*!
  \brief propagates the kinetic polygons in the initialized partition.

  \param k
   maximum number of allowed intersections for each input polygon before its expansion stops.

  \pre successful initialization and k != 0
  */
  void partition(std::size_t k) {

    for (std::size_t idx : m_partitions) {
      Sub_partition& partition = m_partition_nodes[idx];
      Timer timer;
      std::cout.precision(20);

      // Already initialized?
      if (partition.m_data->number_of_support_planes() < 6) {
        std::cout << "Kinetic partition not initialized or empty. Number of support planes: " << partition.m_data->number_of_support_planes() << std::endl;

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

      Propagation propagation(*partition.m_data, m_parameters);
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

      Finalizer finalizer(*partition.m_data, m_parameters);

      if (m_parameters.verbose)
        std::cout << "* getting volumes ..." << std::endl;

      finalizer.create_polyhedra();
      timer.stop();
      const double time_to_finalize = timer.time();

      if (m_parameters.verbose)
        std::cout << "* found all together " << partition.m_data->number_of_volumes() << " volumes" << std::endl;

      if (m_parameters.debug)
        for (std::size_t i = 0; i < partition.m_data->number_of_support_planes(); i++)
          dump_2d_surface_mesh(*partition.m_data, i, partition.m_data->prefix() + "final-surface-mesh-" + std::to_string(i));

      // Timing.
      if (m_parameters.verbose) {
        std::cout << std::endl << "--- TIMING (sec.):" << std::endl;

        std::cout << "* propagation: " << time_to_propagate << std::endl;
        std::cout << "* finalization: " << time_to_finalize << std::endl;
      }
    }

    merge_partitions(0);

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
    return m_data.vertices().size();
  }

  /*!
  \brief returns the number of faces in the kinetic partition.

  \pre successful partition
  */
  std::size_t number_of_faces() const {
    return m_data.face_to_vertices().size();
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
  const std::vector<Point_3>& vertices() const {
    return m_data.vertices();
  }

  /*!
  \brief Vertex indices of face.

  \param face_index
   index of the query face.

  @return
   vector of vertex indices.

  \pre successful partition
  */
  const std::vector<std::size_t>& vertices(std::size_t face_index) const {
    return m_data.face_to_vertices()[face_index];
  }

  /*!
  \brief Face indices of the volume.

  \param volume_index
   index of the query volume.

  @return
   vector of face indices.

  \pre successful partition
  */
  const std::vector<std::size_t>& faces(std::size_t volume_index) const {
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
    return m_data.face_to_volumes()[face_index];
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
    lcc.clear();/*

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
    }*/
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
  void create_bounding_box(
    const FT enlarge_bbox_ratio,
    const bool reorient,
    std::array<Point_3, 8>& bbox) const {

    if (reorient) {
      initialize_optimal_box(bbox);
    }
    else {
      initialize_axis_aligned_box(bbox);
    }

    CGAL_assertion(bbox.size() == 8);

    enlarge_bounding_box(enlarge_bbox_ratio, bbox);

    const auto& minp = bbox.front();
    const auto& maxp = bbox.back();
    if (m_parameters.verbose) {
      std::cout.precision(20);
      std::cout << "* bounding box minp: " << std::fixed <<
        minp.x() << "\t, " << minp.y() << "\t, " << minp.z() << std::endl;
    }
    if (m_parameters.verbose) {
      std::cout.precision(20);
      std::cout << "* bounding box maxp: " << std::fixed <<
        maxp.x() << "\t, " << maxp.y() << "\t, " << maxp.z() << std::endl;
    }
  }

  void initialize_optimal_box(
    std::array<Point_3, 8>& bbox) const {
    /*

    // Number of input points.
    std::size_t num_points = 0;
    for (const auto& poly : m_input_polygons) {
      num_points += poly.size();
    }

    // Set points.
    std::vector<Point_3> points;
    points.reserve(num_points);
    for (const auto& poly : m_input_polygons) {
      for (const auto& point : poly) {
        const Point_3 ipoint(
          static_cast<FT>(CGAL::to_double(point.x())),
          static_cast<FT>(CGAL::to_double(point.y())),
          static_cast<FT>(CGAL::to_double(point.z())));
        points.push_back(ipoint);
      }
    }

    // Compute optimal bbox.
    // The order of faces corresponds to the standard order from here:
    // https://doc.cgal.org/latest/BGL/group__PkgBGLHelperFct.html#gad9df350e98780f0c213046d8a257358e
    const OBB_traits obb_traits;
    std::array<Point_3, 8> ibbox;
    CGAL::oriented_bounding_box(
      points, ibbox,
      CGAL::parameters::use_convex_hull(true).
      geom_traits(obb_traits));

    for (std::size_t i = 0; i < 8; ++i) {
      const auto& ipoint = ibbox[i];
      const Point_3 point(
        static_cast<FT>(ipoint.x()),
        static_cast<FT>(ipoint.y()),
        static_cast<FT>(ipoint.z()));
      bbox[i] = point;
    }

    const FT bbox_length_1 = KSR::distance(bbox[0], bbox[1]);
    const FT bbox_length_2 = KSR::distance(bbox[0], bbox[3]);
    const FT bbox_length_3 = KSR::distance(bbox[0], bbox[5]);
    CGAL_assertion(bbox_length_1 >= FT(0));
    CGAL_assertion(bbox_length_2 >= FT(0));
    CGAL_assertion(bbox_length_3 >= FT(0));
    const FT tol = KSR::tolerance<FT>();
    if (bbox_length_1 < tol || bbox_length_2 < tol || bbox_length_3 < tol) {
      if (m_parameters.verbose) {
        std::cout << "* warning: optimal bounding box is flat, reverting ..." << std::endl;
      }
      initialize_axis_aligned_box(bbox);
    }
    else {
      if (m_parameters.verbose) {
        std::cout << "* using optimal bounding box" << std::endl;
      }
    }*/
  }

  void initialize_axis_aligned_box(
    std::array<Point_3, 8>& bbox) const {


    Bbox_3 box;
    for (const auto& poly : m_input_polygons) {
      box += CGAL::bbox_3(poly.begin(), poly.end());
    }

    // The order of faces corresponds to the standard order from here:
    // https://doc.cgal.org/latest/BGL/group__PkgBGLHelperFct.html#gad9df350e98780f0c213046d8a257358e
    bbox = {
      Point_3(box.xmin(), box.ymin(), box.zmin()),
      Point_3(box.xmax(), box.ymin(), box.zmin()),
      Point_3(box.xmax(), box.ymax(), box.zmin()),
      Point_3(box.xmin(), box.ymax(), box.zmin()),
      Point_3(box.xmin(), box.ymax(), box.zmax()),
      Point_3(box.xmin(), box.ymin(), box.zmax()),
      Point_3(box.xmax(), box.ymin(), box.zmax()),
      Point_3(box.xmax(), box.ymax(), box.zmax()) };

    const FT bbox_length_1 = KSR::distance(bbox[0], bbox[1]);
    const FT bbox_length_2 = KSR::distance(bbox[0], bbox[3]);
    const FT bbox_length_3 = KSR::distance(bbox[0], bbox[5]);
    CGAL_assertion(bbox_length_1 >= FT(0));
    CGAL_assertion(bbox_length_2 >= FT(0));
    CGAL_assertion(bbox_length_3 >= FT(0));
    const FT tol = KSR::tolerance<FT>();
    if (bbox_length_1 < tol || bbox_length_2 < tol || bbox_length_3 < tol) {
      const FT d = 0.1;

      if (bbox_length_1 < tol) { // yz case
        CGAL_assertion_msg(bbox_length_2 >= tol, "ERROR: DEGENERATED INPUT POLYGONS!");
        CGAL_assertion_msg(bbox_length_3 >= tol, "ERROR: DEGENERATED INPUT POLYGONS!");

        bbox[0] = Point_3(bbox[0].x() - d, bbox[0].y() - d, bbox[0].z() - d);
        bbox[3] = Point_3(bbox[3].x() - d, bbox[3].y() + d, bbox[3].z() - d);
        bbox[4] = Point_3(bbox[4].x() - d, bbox[4].y() + d, bbox[4].z() + d);
        bbox[5] = Point_3(bbox[5].x() - d, bbox[5].y() - d, bbox[5].z() + d);

        bbox[1] = Point_3(bbox[1].x() + d, bbox[1].y() - d, bbox[1].z() - d);
        bbox[2] = Point_3(bbox[2].x() + d, bbox[2].y() + d, bbox[2].z() - d);
        bbox[7] = Point_3(bbox[7].x() + d, bbox[7].y() + d, bbox[7].z() + d);
        bbox[6] = Point_3(bbox[6].x() + d, bbox[6].y() - d, bbox[6].z() + d);
        if (m_parameters.verbose) {
          std::cout << "* setting x-based flat axis-aligned bounding box" << std::endl;
        }

      }
      else if (bbox_length_2 < tol) { // xz case
        CGAL_assertion_msg(bbox_length_1 >= tol, "ERROR: DEGENERATED INPUT POLYGONS!");
        CGAL_assertion_msg(bbox_length_3 >= tol, "ERROR: DEGENERATED INPUT POLYGONS!");

        bbox[0] = Point_3(bbox[0].x() - d, bbox[0].y() - d, bbox[0].z() - d);
        bbox[1] = Point_3(bbox[1].x() + d, bbox[1].y() - d, bbox[1].z() - d);
        bbox[6] = Point_3(bbox[6].x() + d, bbox[6].y() - d, bbox[6].z() + d);
        bbox[5] = Point_3(bbox[5].x() - d, bbox[5].y() - d, bbox[5].z() + d);

        bbox[3] = Point_3(bbox[3].x() - d, bbox[3].y() + d, bbox[3].z() - d);
        bbox[2] = Point_3(bbox[2].x() + d, bbox[2].y() + d, bbox[2].z() - d);
        bbox[7] = Point_3(bbox[7].x() + d, bbox[7].y() + d, bbox[7].z() + d);
        bbox[4] = Point_3(bbox[4].x() - d, bbox[4].y() + d, bbox[4].z() + d);
        if (m_parameters.verbose) {
          std::cout << "* setting y-based flat axis-aligned bounding box" << std::endl;
        }

      }
      else if (bbox_length_3 < tol) { // xy case
        CGAL_assertion_msg(bbox_length_1 >= tol, "ERROR: DEGENERATED INPUT POLYGONS!");
        CGAL_assertion_msg(bbox_length_2 >= tol, "ERROR: DEGENERATED INPUT POLYGONS!");

        bbox[0] = Point_3(bbox[0].x() - d, bbox[0].y() - d, bbox[0].z() - d);
        bbox[1] = Point_3(bbox[1].x() + d, bbox[1].y() - d, bbox[1].z() - d);
        bbox[2] = Point_3(bbox[2].x() + d, bbox[2].y() + d, bbox[2].z() - d);
        bbox[3] = Point_3(bbox[3].x() - d, bbox[3].y() + d, bbox[3].z() - d);

        bbox[5] = Point_3(bbox[5].x() - d, bbox[5].y() - d, bbox[5].z() + d);
        bbox[6] = Point_3(bbox[6].x() + d, bbox[6].y() - d, bbox[6].z() + d);
        bbox[7] = Point_3(bbox[7].x() + d, bbox[7].y() + d, bbox[7].z() + d);
        bbox[4] = Point_3(bbox[4].x() - d, bbox[4].y() + d, bbox[4].z() + d);
        if (m_parameters.verbose) {
          std::cout << "* setting z-based flat axis-aligned bounding box" << std::endl;
        }

      }
      else {
        CGAL_assertion_msg(false, "ERROR: WRONG CASE!");
      }
    }
    else {
      if (m_parameters.verbose) {
        std::cout << "* using axis-aligned bounding box" << std::endl;
      }
    }
  }

  void enlarge_bounding_box(
    const FT enlarge_bbox_ratio,
    std::array<Point_3, 8>& bbox) const {

    FT enlarge_ratio = enlarge_bbox_ratio;
    const FT tol = KSR::tolerance<FT>();
    if (enlarge_bbox_ratio == FT(1)) {
      enlarge_ratio += FT(2) * tol;
    }

    const auto a = CGAL::centroid(bbox.begin(), bbox.end());
    Transform_3 scale(CGAL::Scaling(), enlarge_ratio);
    for (auto& point : bbox)
      point = scale.transform(point);

    const auto b = CGAL::centroid(bbox.begin(), bbox.end());
    Transform_3 translate(CGAL::Translation(), a - b);
    for (auto& point : bbox)
      point = translate.transform(point);
  }

  std::pair<int, int> make_canonical_pair(int i, int j)
  {
    if (i > j) return std::make_pair(j, i);
    return std::make_pair(i, j);
  }

  double build_cdt(CDTplus& cdt, const std::vector<typename Intersection_kernel::Point_3>& points, std::vector<std::vector<std::size_t> > &faces, const typename Intersection_kernel::Plane_3& plane) {
    double area = 0;
    From_exact from_exact;
    To_exact to_exact;

    //check orientation of faces so that they are ccw oriented
    for (std::size_t i = 0; i < faces.size(); ++i) {
      auto& v = faces[i];

      std::size_t j = 0;

      CGAL::Orientation res = CGAL::COLLINEAR;
      bool pos = false;
      bool neg = false;

      const std::string vfilename = std::to_string(i) + ".polylines.txt";
      std::ofstream vout(vfilename);
      vout.precision(20);
      vout << std::to_string(v.size() + 1);
      for (auto p : v) {
        vout << " " << from_exact(points[p]);
      }
      vout << " " << from_exact(points[v[0]]);
      vout << std::endl;
      vout.close();

      for (std::size_t j = 0; j < v.size(); j++) {
        std::size_t k = (j + 1) % v.size();
        std::size_t l = (k + 1) % v.size();

        Point_2 pj = from_exact(plane.to_2d(points[v[j]]));
        Point_2 pk = from_exact(plane.to_2d(points[v[k]]));
        Point_2 pl = from_exact(plane.to_2d(points[v[l]]));

        res = orientation(plane.to_2d(points[v[j]]), plane.to_2d(points[v[k]]), plane.to_2d(points[v[l]]));
        if (res == CGAL::LEFT_TURN)
          pos = true;
        if (res == CGAL::RIGHT_TURN)
          neg = true;
      }

      if (pos && neg)
        std::cout << "face is not convex" << std::endl;

      if (!pos && !neg)
        std::cout << "face is degenerated" << std::endl;

      if (neg)
        std::reverse(v.begin(), v.end());
    }

    std::map<std::size_t, std::size_t> face2vtx, vtx2face;
    std::vector<Vertex_handle> vertices;
    for (auto f : faces)
      for (auto v : f) {
        vertices.push_back(cdt.insert(to_exact(from_exact(plane.to_2d(points[v])))));
        vertices.back()->info().set_index(v);
        vertices.back()->info().set_point(points[v]);
        face2vtx[v] = vertices.size() - 1;
        vtx2face[vertices.size() - 1] = v;
      }

    // TODO: insert a range, but keep the vertices in order
    typedef std::set<std::pair<int, int> > Edges;
    Edges edges;
    typedef std::map<std::pair<Vertex_handle, Vertex_handle>, int > HalfEdges;
    HalfEdges halfedges;
    for (std::size_t i = 0; i < faces.size(); ++i) {
      auto& v = faces[i];
      for (std::size_t j = 0; j < v.size(); ++j) {
        int vj = face2vtx[v[j]];
        int vjj = face2vtx[v[(j + 1) % v.size()]];
        std::pair<Edges::iterator, bool> res = edges.insert(make_canonical_pair(vj, vjj));
#ifdef OVERLAY_2_DEBUG
        int vjjj = face2vtx[v[(j + 2) % v.size()]];
        if (orientation(vertices[vj]->point(), vertices[vjj]->point(), vertices[vjjj]->point()) != CGAL::LEFT_TURN) {
          std::cerr << "orientation( " << vertices[vj]->point() << ", " << vertices[vjj]->point() << ", " << vertices[vjjj]->point() << std::endl;
          std::cerr << orientation(vertices[vj]->point(), vertices[vjj]->point(), vertices[vjjj]->point()) << std::endl;
        }
#endif
        halfedges[std::make_pair(vertices[vj], vertices[vjj])] = i;
        if (res.second) {
          cdt.insert_constraint(vertices[vj], vertices[vjj]);
        }
      }
    }
    for (CDTplus::Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end();++fit) {
#ifdef OVERLAY_2_CHECK
      Point_2 p = from_exact(fit->vertex(0)->point());
      Point_2 q = from_exact(fit->vertex(1)->point());
      Point_2 r = from_exact(fit->vertex(2)->point());
      area += CGAL::area(p, q, r);
#endif
      for (int i = 0; i < 3; i++) {
        CDTplus::Edge e(fit, i);
        HalfEdges::iterator it = halfedges.find(std::make_pair(fit->vertex(CDTplus::ccw(i)), fit->vertex(CDTplus::cw(i))));
        if (it != halfedges.end()) {
          fit->info().id = it->second;
        }
      }
    }

    return area;
  }

  std::pair<double, double> overlay(CDTplus& cdtC, const CDTplus& cdtA, const CDTplus& cdtB, const typename Intersection_kernel::Plane_3& plane) {
    From_exact from_exact;
    To_exact to_exact;
    std::pair<double, double> result;
    cdtC = cdtA;
    typename CDTplus::Constraint_iterator ci = cdtB.constraints_begin();

    std::vector<Vertex_handle> vertices;
    vertices.reserve(2);
    //std::size_t idx = 0;
    for (typename CDTplus::Constraint_iterator ci = cdtB.constraints_begin(); ci != cdtB.constraints_end(); ++ci) {
      for (typename CDTplus::Vertices_in_constraint_iterator vi = cdtB.vertices_in_constraint_begin(*ci); vi != cdtB.vertices_in_constraint_end(*ci); vi++) {
        vertices.push_back(*vi);
      }
/*
      const std::string vfilename = std::to_string(idx) + "-constraint.polylines.txt";
      std::ofstream vout(vfilename);
      vout.precision(20);
      vout << 2;
      vout << " " << from_exact(plane.to_3d(vertices[0]->point()));
      vout << " " << from_exact(plane.to_3d(vertices[1]->point()));
      vout << std::endl;
      vout.close();

      idx++;*/

      if (vertices.size() > 2)
        std::cout << "constraint contains more than 2 vertices!" << std::endl;

      for (std::size_t i = 1; i < vertices.size(); i++) {
        cdtC.insert_constraint(to_exact(from_exact(vertices[i - 1]->point())), to_exact(from_exact(vertices[i]->point())));
      }
      vertices.clear();
    }

    // Generate 3D points corresponding to the intersections
    for (typename CDTplus::Finite_vertices_iterator vit = cdtC.finite_vertices_begin(); vit != cdtC.finite_vertices_end(); ++vit) {
      if (!vit->info().input) {
        vit->info().point_3 = plane.to_3d(vit->point());
      }
    }

    // TODO: collect the centroids, perform Hilbert sort and locate
    // with the previous location as hint where to start
    for (typename CDTplus::Finite_faces_iterator cit = cdtC.finite_faces_begin(); cit != cdtC.finite_faces_end(); ++cit) {
      double a = 0;
#ifdef OVERLAY_2_CHECK
      Point_2 ap = from_exact(cit->vertex(0)->point());
      Point_2 aq = from_exact(cit->vertex(1)->point());
      Point_2 ar = from_exact(cit->vertex(2)->point());
      a = CGAL::area(ap, aq, ar);
#endif
      typename Intersection_kernel::Point_2 p = CGAL::centroid(cit->vertex(0)->point(), cit->vertex(1)->point(), cit->vertex(2)->point());
      typename CDTplus::Face_handle fhA = cdtA.locate(p);

      // if(! cdtA.is_infinite(fhA)){
      if (fhA->info().id != -1) {
        cit->info().idA = fhA->info().id;
        // std::cerr << "A: " << fhA->info().id << std::endl;
        result.first += a;
      }
      typename CDTplus::Face_handle fhB = cdtB.locate(p);
      //    if(! cdtB.is_infinite(fhB)){
      if (fhB->info().id != -1) {
        cit->info().idB = fhB->info().id;
        // std::cerr << "B: " << fhB->info().id << std::endl;
        result.second += a;
      }
    }

    return result;
  }

  void collect_faces(std::size_t partition_idx, std::size_t sp_idx, std::vector<std::pair<std::size_t, std::size_t> >& face_idx, std::vector<std::vector<size_t> >& faces) {
    Sub_partition& p = m_partition_nodes[partition_idx];

    for (std::size_t i = 0; i < p.m_data->volumes().size(); i++) {
      typename Data_structure::Volume_cell& v = p.m_data->volumes()[i];
      for (std::size_t j = 0; j < v.faces.size(); j++) {
        if (v.pfaces[j].first == sp_idx) {
          face_idx.push_back(std::make_pair(i, j));
          faces.push_back(p.m_data->face_to_vertices()[v.faces[j]]);
        }
      }
    }
  }

  void check_faces(std::size_t partition_idx, std::size_t sp_idx) {
    Sub_partition& p = m_partition_nodes[partition_idx];

    for (std::size_t i = 0; i < p.m_data->volumes().size(); i++) {
      typename Data_structure::Volume_cell& v = p.m_data->volumes()[i];
      for (std::size_t j = 0; j < v.faces.size(); j++) {
        if (v.pfaces[j].first == sp_idx) {
          if (v.neighbors[j] == -1)
            std::cout << "neighbor not set partition: " << partition_idx << " volume: " << i << " face: " << j << std::endl;
          else
            std::cout << "neighbor is set partition: " << partition_idx << " volume: " << i << " face: " << j << std::endl;
        }
      }
    }
  }

  void collect_planes(std::size_t partition_idx, std::size_t sp_idx, std::vector<std::size_t> &planes) {
    Sub_partition& part = m_partition_nodes[partition_idx];
    typename Data_structure::Support_plane& sp = part.m_data->support_plane(sp_idx);
    auto pedges = sp.mesh().edges();

    std::set<std::size_t> pl;

    for (auto edge : pedges) {
      if (sp.has_iedge(edge)) {
        for (std::size_t p : part.m_data->intersected_planes(sp.iedge(edge)))
          if (!part.m_data->is_bbox_support_plane(p))
            pl.insert(p);
      }
    }

    std::copy(pl.begin(), pl.end(), std::back_inserter(planes));
  }

  void merge_partitions(std::size_t idx) {
    From_exact from_exact;
    if (!m_partition_nodes[idx].children.empty()) {
      std::size_t lower_idx = m_partition_nodes[idx].children[0];
      std::size_t upper_idx = m_partition_nodes[idx].children[1];

      Sub_partition& lower = m_partition_nodes[lower_idx];
      Sub_partition& upper = m_partition_nodes[upper_idx];

      std::size_t lower_sp = lower.split_plane;
      std::size_t upper_sp = upper.split_plane;

      // Collect faces and volumes that are on that plane (use the neighboring index in volumes)
      std::vector<std::pair<std::size_t, std::size_t> > face_idx_lower, face_idx_upper;
      std::vector<std::vector<std::size_t> >  faces_lower, faces_upper;
      std::vector<typename Intersection_kernel::Point_3> vertices_lower = lower.m_data->exact_vertices(), vertices_upper = upper.m_data->exact_vertices();

      collect_faces(lower_idx, lower_sp, face_idx_lower, faces_lower);
      collect_faces(upper_idx, upper_sp, face_idx_upper, faces_upper);

      /*
            const std::string vfilename = "lower.xyz";
            std::ofstream vout(vfilename);
            vout.precision(20);
            for (std::vector<std::size_t>& f : faces_lower)
              for (std::size_t p : f)
              vout << " " << vertices_lower[p] << std::endl;
            vout << std::endl;
            vout.close();


            const std::string vfilename2 = "upper.xyz";
            std::ofstream vout2(vfilename2);
            vout2.precision(20);
            for (auto f : faces_upper)
              for (auto p : f)
                vout2 << " " << vertices_upper[p] << std::endl;

            vout2 << std::endl;
            vout2.close();*/

      typename Intersection_kernel::Plane_3 plane = lower.m_data->support_plane(lower_sp).exact_plane();

      // Collect Plane_3 from support planes (two vectors, one for each other)
      std::vector<std::size_t> planes_lower, planes_upper;
      collect_planes(lower_idx, lower_sp, planes_lower);
      collect_planes(upper_idx, upper_sp, planes_upper);

      // Remove common planes
      auto lower_it = planes_lower.begin();
      auto upper_it = planes_upper.begin();
      while (lower_it != planes_lower.end()) {
        if (*lower_it == *upper_it) {
          lower_it = planes_lower.erase(lower_it);
          upper_it = planes_upper.erase(upper_it);
          if (upper_it == planes_upper.end())
            break;
        }
        else if (*lower_it < *upper_it) {
          lower_it++;
        }
        else if (*lower_it > *upper_it) {
          upper_it++;
          if (upper_it == planes_upper.end())
            break;
        }
      }

      if (!planes_upper.empty()) {
        split_faces(lower_idx, upper_idx, lower_sp, faces_lower, face_idx_lower, vertices_lower, planes_upper);
      }

      if (!planes_lower.empty()) {
        split_faces(upper_idx, lower_idx, upper_sp, faces_upper, face_idx_upper, vertices_upper, planes_lower);
      }

      // How to merge the two Data_structures
      // Identification of common vertices
      // using support planes to check whether I need to check for common vertices? Seems difficult as every vertex is at the intersection of at least three support planes?

      CDTplus lowerCDT, upperCDT;
      double lower_area = build_cdt(lowerCDT, vertices_lower, faces_lower, plane);
      double upper_area = build_cdt(upperCDT, vertices_upper, faces_upper, plane);

      CDTplus combined;
      std::pair<double, double> areas = overlay(combined, lowerCDT, upperCDT, plane);

      for (std::size_t i = 0; i < faces_lower.size(); i++) {
        typename Data_structure::Volume_cell& v = lower.m_data->volumes()[face_idx_lower[i].first];

        std::vector<typename Intersection_kernel::Point_3> pts;
        pts.reserve(faces_lower[i].size());
        for (std::size_t idx : faces_lower[i])
          pts.push_back(vertices_lower[idx]);

        typename Intersection_kernel::Point_3 c = CGAL::centroid(pts.begin(), pts.end(), CGAL::Dimension_tag<0>());

        typename CDTplus::Face_handle neighbor = upperCDT.locate(plane.to_2d(c));
        if (neighbor->info().id < faces_upper.size()) {
          std::cout << "index " << i << " of lower set to " << face_idx_upper[neighbor->info().id].first << std::endl;
          v.neighbors[face_idx_lower[i].second] = face_idx_upper[neighbor->info().id].first;
        }
        else std::cout << "neighbor of face " << i << " of lower has neighbor " << face_idx_upper[neighbor->info().id].first << " in upper" << std::endl;
      }

      check_faces(lower_idx, lower_sp);

      for (std::size_t i = 0; i < faces_upper.size(); i++) {
        typename Data_structure::Volume_cell& v = upper.m_data->volumes()[face_idx_upper[i].first];

        std::vector<typename Intersection_kernel::Point_3> pts;
        pts.reserve(faces_upper[i].size());
        for (std::size_t idx : faces_upper[i])
          pts.push_back(vertices_upper[idx]);

        typename Intersection_kernel::Point_3 c = CGAL::centroid(pts.begin(), pts.end(), CGAL::Dimension_tag<0>());

        typename CDTplus::Face_handle neighbor = lowerCDT.locate(plane.to_2d(c));
        if (neighbor->info().id < faces_upper.size()) {
          std::cout << "index " << i << " of upper set to " << face_idx_lower[neighbor->info().id].first << std::endl;
          v.neighbors[face_idx_upper[i].second] = face_idx_lower[neighbor->info().id].first;
        }
        else std::cout << "neighbor of face " << i << " of upper has neighbor " << face_idx_lower[neighbor->info().id].first << " in upper" << std::endl;
      }

      check_faces(upper_idx, upper_sp);
    }
  }

  void split_faces(std::size_t idx, std::size_t other, std::size_t sp_idx, std::vector<std::vector<std::size_t> > &faces, std::vector<std::pair<std::size_t, std::size_t> >& face_idx, std::vector<typename Intersection_kernel::Point_3> &vertices, std::vector<std::size_t> &planes) {
    typename Intersection_kernel::Plane_3 plane = m_partition_nodes[idx].m_data->support_plane(sp_idx).exact_plane();
    std::vector<typename Intersection_kernel::Point_2> v2d(vertices.size());
    for (std::size_t i = 0; i < vertices.size(); i++)
      v2d[i] = plane.to_2d(vertices[i]);

    From_exact from_exact;

    for (std::size_t pl : planes) {
      typename Intersection_kernel::Line_3 line;
      bool intersect = Data_structure::intersection(plane, m_partition_nodes[other].m_data->support_plane(pl).exact_plane(), line);
      CGAL_assertion(intersect);
      typename Intersection_kernel::Line_2 l2 = m_partition_nodes[idx].m_data->support_plane(sp_idx).to_2d(line);
      //typename Kernel::Line_2 l2 = from_exact(l2_exact);

      std::size_t num_faces = faces.size();

      for (std::size_t f = 0; f < faces.size(); f++) {
        bool neg = false, pos = false;
        for (std::size_t p : faces[f]) {
          CGAL::Oriented_side s = l2.oriented_side(v2d[p]);
          if (s == CGAL::ON_POSITIVE_SIDE) {
            if (neg) {
              CGAL_assertion(f < num_faces);
              split_face(idx, f, faces, face_idx, v2d, plane, vertices, l2);
              break;
            }
            else pos = true;
          }

          if (s == CGAL::ON_NEGATIVE_SIDE) {
            if (pos) {
              CGAL_assertion(f < num_faces);
              split_face(idx, f, faces, face_idx, v2d, plane, vertices, l2);
              break;
            }
            else neg = true;
          }
        }
      }
    }
  }

  void split_face(std::size_t partition, std::size_t f, std::vector<std::vector<std::size_t> >& faces, std::vector<std::pair<std::size_t, std::size_t> > &face_idx, std::vector<typename Intersection_kernel::Point_2> &v2d, typename Intersection_kernel::Plane_3 &plane, std::vector<typename Intersection_kernel::Point_3> &pts, const typename Intersection_kernel::Line_2 &line) {
    std::vector<std::size_t> pos, neg;
    From_exact from_exact;

    const std::string vfilename = std::to_string(f) + "-before.polylines.txt";
    std::ofstream vout(vfilename);
    vout.precision(20);
    vout << std::to_string(faces[f].size() + 1);
    for (auto p : faces[f]) {
      vout << " " << from_exact(pts[p]);
    }
    vout << " " << from_exact(pts[faces[f][0]]);
    vout << std::endl;
    vout.close();

    CGAL::Oriented_side former = line.oriented_side(v2d[faces[f][0]]);

    if (former == CGAL::ON_POSITIVE_SIDE || former== CGAL::ON_ORIENTED_BOUNDARY)
      pos.push_back(faces[f][0]);

    if (former == CGAL::ON_NEGATIVE_SIDE || former == CGAL::ON_ORIENTED_BOUNDARY)
      neg.push_back(faces[f][0]);

    for (std::size_t i = 1; i < faces[f].size() + 1; i++) {
      // Wrap around index
      std::size_t idx = i % faces[f].size();
      CGAL::Oriented_side cur = line.oriented_side(v2d[faces[f][idx]]);
      if (cur == CGAL::ON_ORIENTED_BOUNDARY) {
        neg.push_back(faces[f][idx]);
        pos.push_back(faces[f][idx]);
        former = cur;
        continue;
      }

      // Switching sides without stepping on the line.
      if (cur != former && cur != CGAL::ON_ORIENTED_BOUNDARY && former != CGAL::ON_ORIENTED_BOUNDARY) {
        typename Intersection_kernel::Point_2 p;
        bool intersect = Data_structure::intersection(typename Intersection_kernel::Line_2(v2d[faces[f][idx]], v2d[faces[f][i - 1]]), line, p);
        v2d.push_back(p);
        pts.push_back(plane.to_3d(p));
        pos.push_back(v2d.size() - 1);
        neg.push_back(v2d.size() - 1);
      }

      if (cur == CGAL::ON_POSITIVE_SIDE) {
        pos.push_back(faces[f][idx]);
      }

      if (cur == CGAL::ON_NEGATIVE_SIDE) {
        neg.push_back(faces[f][idx]);
      }

      former = cur;
    }

    bool replaced = false;
    auto& face2vertices = m_partition_nodes[partition].m_data->face_to_vertices();
    auto& volumes = m_partition_nodes[partition].m_data->volumes();

    if (neg.size() > 2) {
      // Check collinearity
      bool collinear = true;
      for (std::size_t i = 2; i < neg.size(); i++) {
        if (!CGAL::collinear(v2d[neg[0]], v2d[neg[1]], v2d[neg[i]])) {
          collinear = false;
          break;
        }
      }

      faces[f] = neg;
      face2vertices[volumes[face_idx[f].first].faces[face_idx[f].second]] = neg;
      replaced = true;
    }

    if (pos.size() > 2) {
      // Check collinearity
      bool collinear = true;
      for (std::size_t i = 2; i < pos.size(); i++) {
        if (!CGAL::collinear(v2d[pos[0]], v2d[pos[1]], v2d[pos[i]])) {
          collinear = false;
          break;
        }
      }
      if (replaced) {
        faces.push_back(pos);
        face2vertices.push_back(pos);
        volumes[face_idx[f].first].faces.push_back(face2vertices.size());
        volumes[face_idx[f].first].neighbors.push_back(volumes[face_idx[f].first].neighbors[face_idx[f].second]);
        volumes[face_idx[f].first].pfaces.push_back(volumes[face_idx[f].first].pfaces[face_idx[f].second]);
        volumes[face_idx[f].first].pface_oriented_outwards.push_back(volumes[face_idx[f].first].pface_oriented_outwards[face_idx[f].second]);
        face_idx.push_back(std::make_pair(face_idx[f].first, volumes[face_idx[f].first].faces.size() - 1));

        m_partition_nodes[partition].m_data->face_to_support_plane().push_back(-1);
        m_partition_nodes[partition].m_data->face_to_volumes().push_back(std::make_pair(-1, -1));
      }
      else {
        faces[f] = pos;
        face2vertices[volumes[face_idx[f].first].faces[face_idx[f].second]] = pos;
      }
    }

    const std::string vfilename2 = std::to_string(f) + "-pos.polylines.txt";
    std::ofstream vout2(vfilename2);
    vout2.precision(20);
    vout2 << std::to_string(pos.size() + 1);
    for (auto p : pos) {
      vout2 << " " << from_exact(pts[p]);
    }
    vout2 << " " << from_exact(pts[pos[0]]);
    vout2 << std::endl;
    vout2.close();

    const std::string vfilename3 = std::to_string(f) + "-neg.polylines.txt";
    std::ofstream vout3(vfilename3);
    vout3.precision(20);
    vout3 << std::to_string(neg.size() + 1);
    for (auto p : neg) {
      vout3 << " " << from_exact(pts[p]);
    }
    vout3 << " " << from_exact(pts[neg[0]]);
    vout3 << std::endl;
    vout3.close();
  }

  void split_partition(std::size_t idx) {
    // Assuming the bbox is axis-aligned

    if (m_partition_nodes[idx].parent != -1) {
      m_partitions.push_back(idx);
      return;
    }

    // Create two children
    m_partition_nodes.resize(m_partition_nodes.size() + 2);

    std::size_t lower_y = m_partition_nodes.size() - 2;
    std::size_t upper_y = lower_y + 1;

    m_partition_nodes[idx].children.push_back(lower_y);
    m_partition_nodes[idx].children.push_back(upper_y);

    m_partition_nodes[lower_y].parent = idx;
    m_partition_nodes[upper_y].parent = idx;

    FT split = (m_partition_nodes[idx].bbox[0].y() + m_partition_nodes[idx].bbox[2].y()) * 0.5;

    // Create bbox and fill in support planes
    //partition2bbox[0] = Bbox_3(bbox.xmin(), bbox.ymin(), bbox.zmin(), bbox.xmax(), split, bbox.zmax());
    //partition2bbox[1] = Bbox_3(bbox.xmin(), split, bbox.zmin(), bbox.xmax(), bbox.ymax(), bbox.zmax());

    // Copy 4 bbox corner points on the lower y side to lower_y partition
    m_partition_nodes[lower_y].bbox[0] = m_partition_nodes[idx].bbox[0];
    m_partition_nodes[lower_y].bbox[1] = m_partition_nodes[idx].bbox[1];
    m_partition_nodes[lower_y].bbox[5] = m_partition_nodes[idx].bbox[5];
    m_partition_nodes[lower_y].bbox[6] = m_partition_nodes[idx].bbox[6];

    // Copy 4 bbox corner points on the upper y side to upper_y partition
    m_partition_nodes[upper_y].bbox[2] = m_partition_nodes[idx].bbox[2];
    m_partition_nodes[upper_y].bbox[3] = m_partition_nodes[idx].bbox[3];
    m_partition_nodes[upper_y].bbox[4] = m_partition_nodes[idx].bbox[4];
    m_partition_nodes[upper_y].bbox[7] = m_partition_nodes[idx].bbox[7];

    // Insert new bbox on split plane
    m_partition_nodes[lower_y].bbox[2] = m_partition_nodes[upper_y].bbox[1] = Point_3(m_partition_nodes[idx].bbox[1].x(), split, m_partition_nodes[idx].bbox[1].z());
    m_partition_nodes[lower_y].bbox[3] = m_partition_nodes[upper_y].bbox[0] = Point_3(m_partition_nodes[idx].bbox[3].x(), split, m_partition_nodes[idx].bbox[3].z());
    m_partition_nodes[lower_y].bbox[4] = m_partition_nodes[upper_y].bbox[5] = Point_3(m_partition_nodes[idx].bbox[4].x(), split, m_partition_nodes[idx].bbox[4].z());
    m_partition_nodes[lower_y].bbox[7] = m_partition_nodes[upper_y].bbox[6] = Point_3(m_partition_nodes[idx].bbox[6].x(), split, m_partition_nodes[idx].bbox[6].z());

    for (std::size_t i = 0; i < m_partition_nodes[idx].input_polygons.size(); i++) {
      bool neg = false, pos = false;
      for (const Point_3& p : m_input_polygons[m_partition_nodes[idx].input_polygons[i]]) {
        if (!neg && p.y() < split) {
          neg = true;
          m_partition_nodes[lower_y].input_polygons.push_back(i);
          if (pos)
            break;
        }
        else if (!pos && p.y() > split) {
          pos = true;
          m_partition_nodes[upper_y].input_polygons.push_back(i);
          if (neg)
            break;
        }
      }
    }

    m_partition_nodes[lower_y].split_plane = 3;
    m_partition_nodes[upper_y].split_plane = 1;

    split_partition(lower_y);
    split_partition(upper_y);
  }

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
