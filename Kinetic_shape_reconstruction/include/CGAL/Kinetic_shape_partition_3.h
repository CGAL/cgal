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

// Internal includes.
#include <CGAL/KSR/utils.h>
#include <CGAL/KSR/debug.h>
#include <CGAL/KSR/parameters.h>

#include <CGAL/KSR_3/Data_structure.h>
#include <CGAL/KSR_3/Initializer.h>
#include <CGAL/KSR_3/FacePropagation.h>
#include <CGAL/KSR_3/Finalizer.h>

#include <CGAL/Octree.h>
#include <CGAL/Orthtree_traits_polygons.h>

//#define OVERLAY_2_DEBUG
#define OVERLAY_2_CHECK

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

  using Index = std::pair<std::size_t, std::size_t>;

private:
  using FT = typename Kernel::FT;
  using Point_2 = typename Kernel::Point_2;
  using Vector_2 = typename Kernel::Vector_2;
  using Plane_3 = typename Kernel::Plane_3;
  using Line_3 = typename Kernel::Line_3;
  using Line_2 = typename Kernel::Line_2;
  using Triangle_2 = typename Kernel::Triangle_2;
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
  using Parameters = KSR::Parameters_3<FT>;

  using Octree = CGAL::Orthtree<CGAL::Orthtree_traits_polygons<Kernel> >;
  using Octree_node = typename Octree::Node_index;

  struct VI
  {
    VI()
      : input(false), idx(-1), idx2(-1, -1), idA2(-1, -1), idB2(-1, -1)
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
    std::size_t idx;  // ivertex?
    std::set<Index> adjacent;
    std::set<Index> ids;
    Index idx2, idA2, idB2;
    bool input;
  };

  // Each face gets as id
  // The overlay face gets also the id from A and B
  struct ID {
    ID()
      : id(-1), idA(-1), idB(-1), id2(std::size_t(-1), std::size_t(-1)), idA2(std::size_t(-1), std::size_t(-1)), idB2(std::size_t(-1), std::size_t(-1))
    {}

    ID(const ID& other)
      :id(other.id), idA(other.idA), idB(other.idB), id2(other.id2), idA2(other.idA2), idB2(other.idB2)
    {}

    ID& operator=(const ID& other)
    {
      id = other.id;
      idA = other.idA;
      idB = other.idB;
      id2 = other.id2;
      idA2 = other.idA2;
      idB2 = other.idB2;
      volA = other.volA;
      volB = other.volB;
      return *this;
    }

    int volA, volB;
    Index id2, idA2, idB2;
    int id, idA, idB;
  };

  typedef CGAL::Triangulation_vertex_base_with_info_2<VI, Intersection_kernel> Vbi;
  typedef CGAL::Triangulation_face_base_with_info_2<ID, Intersection_kernel> Fbi;
  typedef CGAL::Constrained_triangulation_face_base_2<Intersection_kernel, Fbi>    Fb;
  typedef CGAL::Triangulation_data_structure_2<Vbi, Fb>       TDS;
  typedef CGAL::Exact_intersections_tag                     Itag;
  typedef CGAL::Constrained_Delaunay_triangulation_2<Intersection_kernel, TDS, Itag> CDT;
  typedef CGAL::Constrained_triangulation_plus_2<CDT>       CDTplus;
  typedef typename CDTplus::Vertex_handle            Vertex_handle;
  typedef typename CDTplus::Face_handle              Face_handle;
  typedef typename CDTplus::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename CDTplus::Finite_faces_iterator    Finite_faces_iterator;

private:
  struct Sub_partition {
    Sub_partition() : parent(-1) {}
    std::shared_ptr<Data_structure> m_data;
    std::array<typename Intersection_kernel::Point_3, 8> bbox;
    std::vector<typename Intersection_kernel::Plane_3> m_bbox_planes;
    std::vector<std::size_t> input_polygons;
    std::vector<std::vector<Point_3> > clipped_polygons;
    std::vector<typename Intersection_kernel::Plane_3> m_input_planes;
    std::size_t parent;
    std::vector<std::size_t> children;
    std::size_t split_plane;
    std::size_t index;

    std::vector<std::pair<Index, Index> > face_neighbors;
    std::vector<std::vector<Index> > face2vertices;

    std::vector<typename Data_structure::Volume_cell> volumes;
    //std::vector<std::vector<std::size_t> > face2vertices;
    //std::vector<Point_3> exact_vertices;

    typename Octree::Node_index node;
  };

  Parameters m_parameters;
  std::array<Point_3, 8> m_bbox;
  std::vector<Sub_partition> m_partition_nodes; // Tree of partitions.
  std::vector<std::size_t> m_partitions; // Contains the indices of the leaf nodes, the actual partitions to be calculated.
  std::size_t m_num_events;
  std::vector<Point_3> m_points;
  std::vector<std::vector<std::size_t> > m_polygons;
  std::vector<std::vector<Point_3> > m_input_polygons;
  std::vector<typename Intersection_kernel::Plane_3> m_input_planes;
  std::vector<Point_2> m_input_centroids;
  std::vector<std::size_t> m_input2regularized; // Mapping from actual input planes to regularized input planes.
  std::vector<std::vector<std::size_t> > m_regularized2input; // Mapping from partitioning planes to original input polygons.
  std::unique_ptr<Octree> m_octree;
  std::vector<std::size_t> m_node2partition;

  std::vector<Index> m_volumes;
  std::map<Index, std::size_t> m_index2volume;

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
    m_num_events(0),
    m_input2regularized()
  {
    m_parameters.angle_tolerance = parameters::choose_parameter(parameters::get_parameter(np, internal_np::angle_tolerance), 5);
    m_parameters.distance_tolerance = parameters::choose_parameter(parameters::get_parameter(np, internal_np::distance_tolerance), 0.05);
  }

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
    m_num_events(0),
    m_input2regularized()
  {
    m_parameters.angle_tolerance = parameters::choose_parameter(parameters::get_parameter(np, internal_np::angle_tolerance), 5);
    m_parameters.distance_tolerance = parameters::choose_parameter(parameters::get_parameter(np, internal_np::distance_tolerance), 0.05);
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
    To_exact to_exact;
    From_exact from_exact;
    std::size_t offset = m_input2regularized.size();
    for (std::size_t p = 0; p < polygon_range.size();p++) {
      auto& poly = polygon_range[p];

      std::vector<Point_3> pts;
      pts.reserve(poly.size());
      for (auto it : poly)
        pts.push_back(*(input_range.begin() + it));
      Plane_3 pl;
      Point_2 c;
      std::vector<Point_2> ch;
      process_input_polygon(pts, pl, c, ch);
      typename Intersection_kernel::Plane_3 exact_pl = to_exact(pl);

      bool merge = false;
      std::size_t i;
      for (i = 0;i<m_input_planes.size();i++)
        if (within_tolerance(from_exact(m_input_planes[i]), m_input_centroids[i], pl, c)) {
          merge = true;
          break;
        }

      if (merge) {
        m_input2regularized.push_back(i);
        m_regularized2input[i].push_back(p);
        // How to merge? Just do a linear least squares fitting on the full set of points?
        m_input_polygons[i].reserve(m_input_polygons[i].size() + ch.size());

        for (std::size_t j = 0; j < ch.size(); j++)
          m_input_polygons[i].push_back(pl.to_3d(ch[j]));

        process_input_polygon(m_input_polygons[i], pl, m_input_centroids[i], ch);
        m_input_planes[i] = to_exact(pl);

        m_input_polygons[i].resize(ch.size());

        // update centroid of merged plane
        for (std::size_t j = 0; j < ch.size(); j++)
          m_input_polygons[i][j] = pl.to_3d(ch[j]);
      }
      else {
        m_input2regularized.push_back(m_input_planes.size());
        m_regularized2input.push_back(std::vector<std::size_t>());
        m_regularized2input.back().push_back(p);
        m_input_planes.push_back(to_exact(pl));
        m_input_centroids.push_back(c);
        m_input_polygons.push_back(std::vector<Point_3>(ch.size()));

        for (std::size_t i = 0; i < ch.size(); i++)
          m_input_polygons.back()[i] = pl.to_3d(ch[i]);
      }
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
      parameters::get_parameter(np, internal_np::bbox_dilation_ratio), FT(12) / FT(10));
    m_parameters.angle_tolerance = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::angle_tolerance), FT(0) / FT(10));
    m_parameters.distance_tolerance = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::distance_tolerance), FT(0) / FT(10));
    m_parameters.reorient_bbox = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::reorient_bbox), false);


    //CGAL_add_named_parameter(max_octree_depth_t, max_octree_depth, max_octree_depth)
    //CGAL_add_named_parameter(max_octree_node_size_t, max_octree_node_size, max_octree_node_size)

    std::cout.precision(20);
    if (m_input_polygons.size() == 0) {
      std::cout << "Warning: Your input is empty!";
      return;
    }

    std::set<std::size_t> n;
    for (auto p : m_input2regularized)
      n.insert(p);

    assert(m_regularized2input.size() == m_input_polygons.size());
    assert(m_regularized2input.size() == n.size());

    //if (m_parameters.verbose)
      std::cout << m_input2regularized.size() << " input polygons regularized into " << m_input_polygons.size() << " input planes" << std::endl;

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

    if (m_parameters.debug) {
      for (std::size_t i = 0; i < m_input_polygons.size(); i++)
        KSR_3::dump_polygon(m_input_polygons[i], std::to_string(i) + "-input_polygon");
    }

    split_octree();
    m_partitions.resize(m_partition_nodes.size());
    std::iota(m_partitions.begin(), m_partitions.end(), 0);

    for (std::size_t idx : m_partitions) {
      Sub_partition& partition = m_partition_nodes[idx];
      std::cout << idx << ". " << partition.input_polygons.size() << " polygons " << std::flush;
      partition.index = idx;

      partition.m_data = std::make_shared<Data_structure>(m_parameters, std::to_string(idx) + "-");

/*
      std::vector<std::vector<Point_3> > input_polygons(partition.input_polygons.size());
      for (std::size_t i = 0; i < partition.input_polygons.size(); i++)
        input_polygons[i] = m_input_polygons[partition.input_polygons[i]];*/

      Initializer initializer(partition.clipped_polygons, partition.m_input_planes, *partition.m_data, m_parameters);
      initializer.initialize(partition.bbox, partition.input_polygons);
      std::cout << std::endl;
    }

    // Timing.
    if (m_parameters.verbose) {
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
    FT a, b, c;
    partition(k, a, b, c);
  }

#ifndef DOXYGEN_RUNNING
  void partition(std::size_t k, FT &partition_time, FT &finalization_time, FT &conformal_time) {
    m_volumes.clear();
    Timer timer;
    timer.start();
    partition_time = 0;
    finalization_time = 0;
    conformal_time = 0;

    for (std::size_t idx : m_partitions) {
      Sub_partition& partition = m_partition_nodes[idx];
      timer.reset();
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
      std::size_t num_queue_calls = 0;

      Propagation propagation(*partition.m_data, m_parameters);
      std::tie(num_queue_calls, m_num_events) = propagation.propagate(k);

      partition_time += timer.time();

      if (m_parameters.verbose) {
        std::cout << "* propagation finished" << std::endl;
        std::cout << "* number of queue calls: " << num_queue_calls << std::endl;
        std::cout << "* number of events handled: " << m_num_events << std::endl;
      }

      if (m_parameters.verbose) {
        std::cout << std::endl << "--- FINALIZING PARTITION:" << std::endl;
      }

      // Finalization.

      if (m_parameters.debug)
        for (std::size_t i = 0; i < partition.m_data->number_of_support_planes(); i++)
          if (!partition.m_data->support_plane(i).mesh().is_valid(true))
            std::cout << i << ". support has an invalid mesh!" << std::endl;

      for (std::size_t i = 6; i < partition.m_data->number_of_support_planes(); i++) {
        bool initial = false;
        typename Data_structure::Support_plane& sp = partition.m_data->support_plane(i);

        for (const auto &f : sp.mesh().faces())
          if (sp.is_initial(f)) {
            initial = true;
            break;
          }

        if (!initial)
          std::cout << i << " sp has no initial face before" << std::endl;
      }

      Finalizer finalizer(*partition.m_data, m_parameters);

      if (m_parameters.verbose)
        std::cout << "* getting volumes ..." << std::endl;

      finalizer.create_polyhedra();
      finalization_time += timer.time();

      for (std::size_t i = 6; i < partition.m_data->number_of_support_planes(); i++) {
        bool initial = false;
        typename Data_structure::Support_plane& sp = partition.m_data->support_plane(i);

        for (const auto& f : sp.mesh().faces())
          if (sp.is_initial(f)) {
            initial = true;
            break;
          }

        if (!initial)
          std::cout << i << " sp has no initial face" << std::endl;
      }

      if (m_parameters.verbose)
        std::cout << idx << ". partition with " << partition.input_polygons.size() << " input polygons split into " << partition.m_data->number_of_volumes() << " volumes" << std::endl;

/*
      if (m_parameters.debug)
        for (std::size_t i = 0; i < partition.m_data->number_of_support_planes(); i++)
          dump_2d_surface_mesh(*partition.m_data, i, partition.m_data->prefix() + "final-surface-mesh-" + std::to_string(i));*/
    }

    //for (std::size_t i = 0;i<m_partition_nodes.size();i++)
    //  std::cout << i << ". partition with " << m_partition_nodes[i].input_polygons.size() << " input polygons split into " << m_partition_nodes[i].m_data->number_of_volumes() << " volumes" << std::endl;

    // Convert face_neighbors to pair<Index, Index>
    for (std::size_t i = 0; i < m_partitions.size(); i++) {
      Sub_partition& partition = m_partition_nodes[m_partitions[i]];

      for (std::size_t j = 0; j < partition.m_data->number_of_volumes(); j++) {
        m_volumes.push_back(std::make_pair(m_partitions[i], j));
      }

      partition.face_neighbors.resize(partition.m_data->face_to_volumes().size());
      for (std::size_t j = 0; j < partition.m_data->face_to_volumes().size(); j++) {
        auto& p = partition.m_data->face_to_volumes()[j];
        partition.face_neighbors[j] = std::make_pair(Index(m_partitions[i], p.first), Index(m_partitions[i], p.second));
      }

      partition.face2vertices.resize(partition.m_data->face_to_vertices().size());

      for (std::size_t j = 0; j < partition.m_data->face_to_vertices().size(); j++) {
        partition.face2vertices[j].resize(partition.m_data->face_to_vertices()[j].size());
        for (std::size_t k = 0; k < partition.m_data->face_to_vertices()[j].size(); k++)
          partition.face2vertices[j][k] = std::make_pair(m_partitions[i], partition.m_data->face_to_vertices()[j][k]);
      }
    }

    for (std::size_t i = 0; i < m_volumes.size(); i++)
      m_index2volume[m_volumes[i]] = i;

    timer.stop();
    timer.reset();
    timer.start();
    make_conformal(0);
    conformal_time = timer.time();

    return;
  }
#endif

  /// @}

  /*******************************
  **         Access         **
  ********************************/

  /// \name Access
  /// @{
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
    return m_volumes.size();
  }

  /*!
  \brief returns barycenter for a given volume.

  \param volume_index
  index of the volume.

  \pre successful partition
  */
  const Point_3 &volume_centroid(std::size_t volume_index) const {
    assert(volume_index < m_volumes.size());
    auto p = m_volumes[volume_index];
    return m_partition_nodes[p.first].m_data->volumes()[p.second].centroid;
  }

  /*!
  \brief provides faces of the partition belonging to the regularized input polygon.

  \tparam OutputIterator
  must be an output iterator to which `Index` can be assigned.

  \param polygon_index
  index of the regularized input polygon.

  \param it
  output iterator.

  \pre successful partition
  */
  template<class OutputIterator>
  void faces_of_regularized_polygon(const std::size_t polygon_index, OutputIterator it) const {
    if (polygon_index >= m_input_planes.size()) {
      assert(false);
    }

    for (std::size_t idx : m_partitions) {
      const Sub_partition& p = m_partition_nodes[idx];
      // Check if it contains this input polygon and get support plane index
      int sp_idx = -1;
      for (std::size_t i = 0; i < p.input_polygons.size(); i++) {
        if (p.input_polygons[i] == polygon_index) {
          sp_idx = p.m_data->support_plane_index(i);
          break;
        }
      }

      // Continue if the partition does not contain this input polygon.
      if (sp_idx == -1)
        continue;

      auto pfaces = p.m_data->pfaces(sp_idx);
      auto f2i = p.m_data->face_to_index();

      for (const auto& f : pfaces) {
        assert(f.first == sp_idx);
        auto fit = f2i.find(f);
        assert(fit != f2i.end());

        *it++ = std::make_pair(idx, fit->second);
      }
    }
  }

  /*!
  \brief maps points onto the faces of the input polygon 'polygon_index' in the partition

  \param polygon_index
  index of the regularized input polygon.

  \param pts
  points to be mapped onto the faces of the partition.

  \param mapping
  resulting mapping vector containing one pair for each face in the partition containing points from pts.

  \pre successful partition
  */
  void map_points_to_regularized_polygons(const std::size_t polygon_index, const std::vector<Point_3>& pts, std::vector<std::pair<Index, std::vector<std::size_t> > > &mapping) {
    std::vector<Index> faces;

    if (polygon_index >= m_input_planes.size()) {
      assert(false);
    }

    for (std::size_t idx : m_partitions) {
      const Sub_partition& p = m_partition_nodes[idx];
      // Check if it contains this input polygon and get support plane index
      std::size_t sp_idx = -1;
      for (std::size_t i = 0; i < p.input_polygons.size(); i++) {
        if (p.input_polygons[i] == polygon_index) {
          sp_idx = p.m_data->support_plane_index(i);
          break;
        }
      }

      // Continue if the partition does not contain this input polygon.
      if (sp_idx == -1)
        continue;

      // Filter points
      From_exact from_exact;
      std::array<FT, 3> bbmin = { from_exact(p.bbox[0][0]), from_exact(p.bbox[0][1]), from_exact(p.bbox[0][2]) };
      std::array<FT, 3> bbmax = { from_exact(p.bbox[7][0]), from_exact(p.bbox[7][1]), from_exact(p.bbox[7][2]) };
      assert(bbmin[0] < bbmax[0]);
      assert(bbmin[1] < bbmax[1]);
      assert(bbmin[2] < bbmax[2]);

      std::vector<Point_2> pts2d;
      std::vector<std::size_t> idx2d;
      auto sp = p.m_data->support_plane(sp_idx);

      for (std::size_t i = 0; i < pts.size(); i++) {
        if (bbmin[0] <= pts[i][0] && pts[i][0] <= bbmax[0]
          && bbmin[1] <= pts[i][1] && pts[i][1] <= bbmax[1]
          && bbmin[2] <= pts[i][2] && pts[i][2] <= bbmax[2]) {
          pts2d.push_back(sp.to_2d(pts[i]));
          idx2d.push_back(i);
        }
      }

      /*auto writer = pts.end()--;
      auto reader = pts.begin();
      while (reader < writer) {
        while ((*reader[0] < bbmin[0] || bbmax[0] < *reader[0]
          || *reader[1] < bbmin[1] || bbmax[1] < *reader[1]
          || *reader[2] < bbmin[2] || bbmax[2] < *reader[2]) && reader < writer)
          reader++;

        while ((*bbmin[0] <= *writer[0] && *writer[0] <= bbmax[0]
          && *bbmin[1] <= *writer[1] && *writer[0] <= bbmax[1]
          && *bbmin[2] <= *writer[2] && *writer[0] <= bbmax[2]) && reader < writer)
          writer--;

        if (writer >= reader)
          break;

        auto tmp = *writer;
        *writer = *reader;
        *reader = tmp;

        reader++;
        writer--;
      };*/

      const auto& initial = p.m_data->face_is_part_of_input_polygon();
      for (std::size_t f = 0; f < p.m_data->face_to_support_plane().size();f++) {
        if (p.m_data->face_to_support_plane()[f] != sp_idx || !initial[f])
          continue;

        mapping.resize(mapping.size() + 1);
        auto& m = mapping.back();
        m.first = Index(idx, f);

        std::vector<Point_3> vts;
        std::vector<Point_2> vts2d;

        vertices(m.first, std::back_inserter(vts));
        vts2d.reserve(vts.size());

        for (const auto& v : vts)
          vts2d.push_back(sp.to_2d(v));

        // Todo: Remove check if vts are ccw
        Polygon_2<Kernel> poly(vts2d.begin(), vts2d.end());
        if (poly.is_clockwise_oriented())
          std::reverse(vts2d.begin(), vts2d.end());

        for (std::size_t i = 0;i<pts2d.size();i++) {
          const auto& pt = pts2d[i];
          bool outside = false;

          // poly, vertices and edges in IFace are oriented ccw
          std::size_t idx = 0;
          for (std::size_t i = 0; i < vts2d.size(); i++) {
            Vector_2 ts = (vts2d[(i + vts2d.size() - 1) % vts2d.size()]) - pt;
            Vector_2 tt = (vts2d[i]) - pt;

            bool ccw = (tt.x() * ts.y() - tt.y() * ts.x()) <= 0;
            if (!ccw) {
              outside = true;
              break;
            }
          }

          if (outside)
            continue;

          m.second.push_back(idx2d[i]);
        }
      }

      // Order of the vertices should be ccw
      /*IFace face = IFace(-1);
      for (auto& f : initial_faces) {
        Face_property& fp = m_data.igraph().face(f);

        for (const auto& p : pts) {

          typename Intersection_kernel::Point_2& p = to_exact(sp.data().centroid);
          bool outside = false;

          // poly, vertices and edges in IFace are oriented ccw
          std::size_t idx = 0;
          for (std::size_t i = 0; i < fp.pts.size(); i++) {
            typename Intersection_kernel::Vector_2 ts = fp.pts[(i + fp.pts.size() - 1) % fp.pts.size()] - p;
            typename Intersection_kernel::Vector_2 tt = fp.pts[i] - p;

            bool ccw = (tt.x() * ts.y() - tt.y() * ts.x()) <= 0;
            if (!ccw) {
              outside = true;
              break;
            }
          }
          if (!outside) {
            if (face == -1)
              face = f;
            else {
              std::cout << "Two faces found for " << sp_idx << " sp, f1 " << face << " f2 " << f << std::endl;
            }
          }
        }
      }
      if (face != -1) {

        if (!m_data.igraph().face(face).part_of_partition) {
          m_data.add_iface_to_mesh(sp_idx, face);
          sp.data().initial_ifaces.push_back(face);
        }
      }
      else
        std::cout << "No IFace found for sp " << sp_idx << std::endl;*/


      /*auto pfaces = p.m_data->pfaces(sp_idx);
      auto f2i = p.m_data->face_to_index();

      for (const auto& f : pfaces) {
        assert(f.first == sp_idx);
        auto fit = f2i.find(f);
        assert(fit != f2i.end());

        *it++ = std::make_pair(idx, fit->second);
      }*/
    }
  }

  /*!
  \brief provides the exact 'Plane_3' for a regularized input polygon.

  \param polygon_index
  index of regularized input polygon.

  */
  const typename Intersection_kernel::Plane_3 &regularized_plane(std::size_t polygon_index) const {
    return m_input_planes[polygon_index];
  }

  /*!
  \brief provides the mapping of regularized input planes to inserted input planes.

  \return a vector containing the indices of input polygons for each regularized input polygon.

  */
  const std::vector<std::vector<std::size_t> > &regularized_input_mapping() const {
    return m_regularized2input;
  }

  /*!
  \brief Mapping of a vertex index to its position.

  \param vertex_index
  query vertex.

  \return 'GeomTraits::Point_3' of a vertex.

  \pre successful partition
  */
  const Point_3& vertex(const Index& vertex_index) const {
    return m_partition_nodes[vertex_index.first].m_data->vertices()[vertex_index.second];
  }

  /*!
  \brief Mapping of a vertex index to its exact position.

  \param vertex_index
  query vertex.

  \return 'Intersection_kernel::Point_3' of a vertex.

  \pre successful partition
  */
  const typename Intersection_kernel::Point_3& exact_vertex(const Index& vertex_index) const {
    return m_partition_nodes[vertex_index.first].m_data->exact_vertices()[vertex_index.second];
  }

  /*!
  \brief Vertex positions of a face of the kinetic partition.

  \tparam OutputIterator
  must be an output iterator to which `GeomTraits::Point_3` can be assigned.

  \param face_index
   index of the query face.

  \param it
  output iterator.

  \pre successful partition
  */
  template<class OutputIterator>
  void vertices(const Index& face_index, OutputIterator it) const {
    for (auto& p : m_partition_nodes[face_index.first].face2vertices[face_index.second])
      *it++ = m_partition_nodes[p.first].m_data->vertices()[p.second];
  }

  /*!
  \brief Vertices of a face of the kinetic partition.

  \tparam OutputIterator
  must be an output iterator to which `Index` can be assigned.

  \param face_index
   index of the query face.

  \param it
  output iterator.

  \pre successful partition
  */
  template<class OutputIterator>
  void vertex_indices(const Index& face_index, OutputIterator it) const {
  for (auto& i : m_partition_nodes[face_index.first].m_data->face_to_vertices()[face_index.second])
    *it++ = std::make_pair(face_index.first, i);
  }

  /*!
  \brief Exact vertex positions of a face of the kinetic partition.

  \tparam OutputIterator
  must be an output iterator to which `Intersection_kernel::Point_3` can be assigned.

  \param face_index
   index of the query face.

  \param it
  output iterator.

  \pre successful partition
  */
  template<class OutputIterator>
  void exact_vertices(const Index& face_index, OutputIterator it) const {

    for (auto& p : m_partition_nodes[face_index.first].face2vertices[face_index.second])
      *it++ = m_partition_nodes[p.first].m_data->exact_vertices()[p.second];
  }

  /*!
  \brief Vertices and their exact positions of a face of the kinetic partition.

  \tparam OutputIterator
  must be an output iterator to which `Intersection_kernel::Point_3` can be assigned.

  \tparam IndexOutputIterator
  must be an output iterator to which `Index` can be assigned.

  \param face_index
   index of the query face.

  \param pit
  output iterator for vertex positions.

  \param iit
  output iterator for vertices.

  \pre successful partition
  */
  template<class OutputIterator, class IndexOutputIterator>
  void exact_vertices(const Index& face_index, OutputIterator pit, IndexOutputIterator iit) const {
    const auto& v = m_partition_nodes[face_index.first].m_data->exact_vertices();
    for (auto& i : m_partition_nodes[face_index.first].m_data->face_to_vertices()[face_index.second]) {
      *iit++ = std::make_pair(face_index.first, i);
      *pit++ = v[i];
    }
  }

  /*!
  \brief Face indices of a volume.

  \tparam OutputIterator
  must be an output iterator to which `Index` can be assigned.

  \param volume_index
   index of the query volume.

  \param it
  output iterator.

  \pre successful partition
  */
  template<class OutputIterator>
  void faces(std::size_t volume_index, OutputIterator it) const {
    CGAL_assertion(m_volumes.size() > volume_index);
    auto p = m_volumes[volume_index];

    for (std::size_t i : m_partition_nodes[p.first].m_data->volumes()[p.second].faces)
      *it++ = std::make_pair(p.first, i);
  }

  /*!
  \brief Retrieves all unique faces of the partition. Unique means that the shared face between two adjacent volumes is only contained once.

  \tparam OutputIterator
  must be an output iterator to which `Index` can be assigned.

  \param it
  output iterator.

  \pre successful partition
  */
  template<class OutputIterator>
  void unique_faces(OutputIterator it) const {
    for (std::size_t i = 0; i < m_partition_nodes.size(); i++) {
      const Sub_partition& p = m_partition_nodes[i];
      for (std::size_t j = 0; j < p.face_neighbors.size(); j++) {
        if (p.face_neighbors[j].second.first == i)
          *it++ = Index(i, j);
        else if (i < p.face_neighbors[j].second.first)
          *it++ = Index(i, j);
      }
    }
  }

  /*!
  \brief Indices of adjacent volumes. Negative indices correspond to the empty spaces behind the sides of the bounding box.

    -1 zmin
    -2 ymin
    -3 xmax
    -4 ymax
    -5 xmin
    -6 zmax

    \param face_index
    index of the query face.

    @return
    pair of adjacent volumes.

    \pre successful partition
  */
  const std::pair<int, int> neighbors(const Index &face_index) const {
    const auto &p = m_partition_nodes[face_index.first].face_neighbors[face_index.second];
    if (p.second.second >= std::size_t(-6)) { // Faces on the boundary box are neighbors with an infinite outside volume
      auto it = m_index2volume.find(p.first);
      assert(it != m_index2volume.end());
      return std::pair<int, int>(static_cast<int>(it->second), static_cast<int>(p.second.second));
    }
    else {
      auto it1 = m_index2volume.find(p.first);
      assert(it1 != m_index2volume.end());
      auto it2 = m_index2volume.find(p.second);
      assert(it2 != m_index2volume.end());
      return std::pair<int, int>(static_cast<int>(it1->second), static_cast<int>(it2->second));
    }
    //const auto& p = m_partition_nodes[face_index.first].m_data->face_to_volumes()[face_index.second];
    //return std::pair<Index, Index>(std::make_pair(face_index.first, p.first), std::make_pair(face_index.first, p.second));// m_data.face_to_volumes()[face_index];
  }

  /*
  \brief Retrieves the support plane generated from regularized input polygon.

  \param input_polygon_index
   index of the input polygon.

  @return
   index into polygon_map provided on initialization.

  \pre successful partition

  std::size_t support_plane_index(const std::size_t input_polygon_index) const {
      const int support_plane_idx = m_data.support_plane_index(input_polygon_index);
      CGAL_assertion(support_plane_idx >= 6);
      return support_plane_idx;
  }*/

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

  void process_input_polygon(const std::vector<Point_3> poly, Plane_3& pl, Point_2& c, std::vector<Point_2>& ch) const {
    CGAL::linear_least_squares_fitting_3(poly.begin(), poly.end(), pl, CGAL::Dimension_tag<0>());

    std::vector<Point_2> pts2d(poly.size());
    for (std::size_t i = 0; i < poly.size(); i++)
      pts2d[i] = pl.to_2d(poly[i]);

    ch.clear();
    CGAL::convex_hull_2(pts2d.begin(), pts2d.end(), std::back_inserter(ch));

    // Centroid
    FT x = 0, y = 0, w = 0;
    for (std::size_t i = 2; i < ch.size(); i++) {
      Triangle_2 tri(ch[0], ch[i - 1], ch[i]);
      w += CGAL::area(ch[0], ch[i - 1], ch[i]);
      Point_2 c = CGAL::centroid(ch[0], ch[i - 1], ch[i]);
      x += c.x() * w;
      y += c.y() * w;
    }

    c = Point_2(x / w, y / w);
  }

  std::pair<int, int> make_canonical_pair(int i, int j)
  {
    if (i > j) return std::make_pair(j, i);
    return std::make_pair(i, j);
  }

  double build_cdt(CDTplus& cdt, std::vector<Index>& faces, const typename Intersection_kernel::Plane_3& plane) {
    double area = 0;
    From_exact from_exact;
    //To_exact to_exact;

    cdt.clear();

    //check orientation of faces so that they are ccw oriented
    std::vector<std::vector<Index> > pts_idx(faces.size());
    std::vector<std::vector<typename Intersection_kernel::Point_3> > pts(faces.size());
    for (std::size_t i = 0; i < faces.size(); ++i) {
      exact_vertices(faces[i], std::back_inserter(pts[i]), std::back_inserter(pts_idx[i]));
      //auto& v = faces[i];

      std::size_t j = 0;

      CGAL::Orientation res = CGAL::COLLINEAR;
      bool pos = false;
      bool neg = false;

      for (std::size_t j = 0; j < pts[i].size(); j++) {
        std::size_t k = (j + 1) % pts[i].size();
        std::size_t l = (k + 1) % pts[i].size();

        Point_2 pj = from_exact(plane.to_2d(pts[i][j]));
        Point_2 pk = from_exact(plane.to_2d(pts[i][k]));
        Point_2 pl = from_exact(plane.to_2d(pts[i][l]));

        res = orientation(plane.to_2d(pts[i][j]), plane.to_2d(pts[i][k]), plane.to_2d(pts[i][l]));
        if (res == CGAL::LEFT_TURN)
          pos = true;
        if (res == CGAL::RIGHT_TURN)
          neg = true;
      }

      if (pos && neg) {
        std::cout << "face is not convex" << std::endl;
        exit(1);
      }

      if (!pos && !neg) {
        std::cout << "face is degenerated" << std::endl;
        exit(1);
      }

      if (neg) {
        std::reverse(pts[i].begin(), pts[i].end());
        std::reverse(pts_idx[i].begin(), pts_idx[i].end());
      }
    }

    std::map<Index, std::size_t> face2vtx;
    std::map<std::size_t, Index> vtx2face;
    std::vector<Vertex_handle> vertices;
    for (std::size_t f = 0; f < faces.size(); f++)
      for (std::size_t v = 0; v < pts_idx[f].size(); v++) {
        //vertices.push_back(cdt.insert(to_exact(from_exact(plane.to_2d(pts[f][v])))));
        vertices.push_back(cdt.insert(plane.to_2d(pts[f][v])));
        vertices.back()->info().idA2 = pts_idx[f][v];
        assert(pts_idx[f][v].first != -1);
        assert(pts_idx[f][v].second != -1);
        vertices.back()->info().adjacent.insert(faces[f]);
        vertices.back()->info().set_point(pts[f][v]);
        face2vtx[pts_idx[f][v]] = vertices.size() - 1;
        vtx2face[vertices.size() - 1] = pts_idx[f][v];
      }

    typedef std::set<std::pair<int, int> > Edges;
    Edges edges;

    for (std::size_t i = 0; i < pts_idx.size(); ++i) {
      auto& v = pts_idx[i];
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
        if (res.second) {
          cdt.insert_constraint(vertices[vj], vertices[vjj]);
        }
      }
    }

    for (CDTplus::Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
#ifdef OVERLAY_2_CHECK
      Point_2 p = from_exact(fit->vertex(0)->point());
      Point_2 q = from_exact(fit->vertex(1)->point());
      Point_2 r = from_exact(fit->vertex(2)->point());
      area += CGAL::area(p, q, r);
#endif

      std::set<Index>& a(fit->vertex(0)->info().adjacent), & b(fit->vertex(1)->info().adjacent), & c(fit->vertex(2)->info().adjacent);

      std::set<Index> res, res2;
      Index common(std::size_t(-1), std::size_t(-1));
      std::set_intersection(a.begin(), a.end(), b.begin(), b.end(), std::inserter(res, res.begin()));
      std::set_intersection(res.begin(), res.end(), c.begin(), c.end(), std::inserter(res2, res2.begin()));

      if (res2.size() != 1) {
        std::cout << "build_cdt: face assignment not unique!" << std::endl;
        const std::string vfilename = "no-face.polylines.txt";
        std::ofstream vout(vfilename);
        vout.precision(20);
        vout << 4;
        vout << " " << from_exact(plane.to_3d(fit->vertex(0)->point()));
        vout << " " << from_exact(plane.to_3d(fit->vertex(1)->point()));
        vout << " " << from_exact(plane.to_3d(fit->vertex(2)->point()));
        vout << " " << from_exact(plane.to_3d(fit->vertex(0)->point()));
        vout << std::endl;
        vout.close();
      }
      else fit->info().id2 = *res2.begin();
    }

    return area;
  }

  double build_cdt(CDTplus& cdt, std::vector<CDTplus>& partitions, const typename Intersection_kernel::Plane_3& plane) {
    if (partitions.size() == 0)
      return 0;

    double area = 0;

    From_exact from_exact;
    //To_exact to_exact;
    cdt = partitions[0];

    for (std::size_t i = 1; i < partitions.size(); i++) {
      std::vector<Vertex_handle> vertices;
      vertices.reserve(6);

      for (typename CDTplus::Constraint_iterator ci = partitions[i].constraints_begin(); ci != partitions[i].constraints_end(); ++ci) {
        for (typename CDTplus::Vertices_in_constraint_iterator vi = partitions[i].vertices_in_constraint_begin(*ci); vi != partitions[i].vertices_in_constraint_end(*ci); vi++) {
          vertices.push_back(*vi);
        }

        // Insert constraints and replacing vertex handles in vector while copying data.
        for (std::size_t i = 0;i<vertices.size();i++) {
          VI tmp = vertices[i]->info();
          vertices[i] = cdt.insert(vertices[i]->point());
          vertices[i]->info() = tmp;
        }

        for (std::size_t i = 1; i < vertices.size(); i++)
          cdt.insert_constraint(vertices[i - 1], vertices[i]);

        vertices.clear();
      }
    }

    // Generate 3D points corresponding to the intersections
    std::size_t newpts = 0;
    for (typename CDTplus::Finite_vertices_iterator vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit) {
      if (!vit->info().input) {
        vit->info().point_3 = plane.to_3d(vit->point());
        vit->info().idA2 = vit->info().idB2 = vit->info().idx2 = Index(-1, -1);
        newpts++;
      }
    }

    //std::cout << newpts << " new vertices added in build_cdt from cdts" << std::endl;

    for (CDTplus::Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
#ifdef OVERLAY_2_CHECK
      Point_2 p = from_exact(fit->vertex(0)->point());
      Point_2 q = from_exact(fit->vertex(1)->point());
      Point_2 r = from_exact(fit->vertex(2)->point());
      area += CGAL::area(p, q, r);
#endif

      Index idx(std::size_t(-1), std::size_t(-1));

      typename Intersection_kernel::Point_2 pt = CGAL::centroid(fit->vertex(0)->point(), fit->vertex(1)->point(), fit->vertex(2)->point());
      for (std::size_t i = 0; i < partitions.size(); i++) {
        typename CDTplus::Face_handle fh = partitions[i].locate(pt);

        if (!partitions[i].is_infinite(fh)) {
          if (fh->info().id2 != std::make_pair(std::size_t(-1), std::size_t(-1)))
            idx = fit->info().id2 = fh->info().id2;
          else
            std::cout << "Face id is missing " << std::endl;
        }
      }

      if (fit->info().id2.first == std::size_t(-1))
        std::cout << "cdt fusion: no id found" << std::endl;
    }

    return area;
  }

  /*

  double build_cdt(CDTplus& cdt, const std::vector<typename Intersection_kernel::Point_3>& points, const std::vector<std::size_t> &volumes, std::vector<std::vector<std::size_t> >& faces, const typename Intersection_kernel::Plane_3& plane) {
    double area = 0;
    From_exact from_exact;
    To_exact to_exact;

    cdt.clear();

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

      if (pos && neg) {
        std::cout << "face is not convex" << std::endl;
        exit(1);
      }

      if (!pos && !neg) {
        std::cout << "face is degenerated" << std::endl;
        exit(1);
      }

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
    for (CDTplus::Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
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
*/

  std::pair<double, double> overlay(CDTplus& cdtC, const CDTplus& cdtA, const CDTplus& cdtB, const typename Intersection_kernel::Plane_3& plane) {
    From_exact from_exact;
    //To_exact to_exact;
    std::pair<double, double> result;
    cdtC = cdtA;

    std::vector<Vertex_handle> vertices;
    vertices.reserve(2);
    //std::size_t idx = 0;
    for (typename CDTplus::Constraint_iterator ci = cdtB.constraints_begin(); ci != cdtB.constraints_end(); ++ci) {
      for (typename CDTplus::Vertices_in_constraint_iterator vi = cdtB.vertices_in_constraint_begin(*ci); vi != cdtB.vertices_in_constraint_end(*ci); vi++) {
        vertices.push_back(*vi);
      }

//       if (vertices.size() > 2) {
//         const std::string vfilename = std::to_string(idx) + "-constraint.polylines.txt";
//         std::ofstream vout(vfilename);
//         vout.precision(20);
//         vout << vertices.size();
//         for (std::size_t i = 0;i<vertices.size();i++)
//           vout << " " << from_exact(plane.to_3d(vertices[i]->point()));
//         vout << std::endl;
//         vout.close();
//       }

      //idx++;

//       if (vertices.size() > 2
//         std::cout << "constraint contains more than 2 vertices!" << std::endl;

        // Insert constraints and replacing vertex handles in vector while copying data.
      for (std::size_t i = 0;i<vertices.size();i++) {
        VI tmp = vertices[i]->info();
        vertices[i] = cdtC.insert(vertices[i]->point());
        vertices[i]->info().idB2 = tmp.idA2;
      }

      for (std::size_t i = 1; i < vertices.size(); i++) {
        cdtC.insert_constraint(((vertices[i - 1])), ((vertices[i])));
      }
      vertices.clear();
    }

    std::size_t newpts = 0;
    // Generate 3D points corresponding to the intersections
    //std::ofstream vout3("newpts.xyz");
    //vout3.precision(20);
    for (typename CDTplus::Finite_vertices_iterator vit = cdtC.finite_vertices_begin(); vit != cdtC.finite_vertices_end(); ++vit) {
      if (!vit->info().input) {
        vit->info().point_3 = plane.to_3d(vit->point());
        vit->info().idA2 = vit->info().idB2 = vit->info().idx2 = Index(-1, -1);
        //vout3 << " " << from_exact(vit->info().point_3) << std::endl;
        newpts++;
      }
    }
    //vout3 << std::endl;
    //vout3.close();

    //std::cout << newpts << " new vertices added in cdt" << std::endl;

/*
    const std::string vfilename = "location_failures.xyz";
    std::ofstream vout(vfilename);
    vout.precision(20);*/

    // TODO: collect the centroids, perform Hilbert sort and locate
    // with the previous location as hint where to start
    for (typename CDTplus::Finite_faces_iterator cit = cdtC.finite_faces_begin(); cit != cdtC.finite_faces_end(); ++cit) {
      double a = 0;
      cit->info().id2 = std::make_pair(-1, -1);
#ifdef OVERLAY_2_CHECK
      Point_2 ap = from_exact(cit->vertex(0)->point());
      Point_2 aq = from_exact(cit->vertex(1)->point());
      Point_2 ar = from_exact(cit->vertex(2)->point());
      a = CGAL::area(ap, aq, ar);
#endif
      typename Intersection_kernel::Point_2 p = CGAL::centroid(cit->vertex(0)->point(), cit->vertex(1)->point(), cit->vertex(2)->point());
      typename CDTplus::Face_handle fhA = cdtA.locate(p);

      if (cdtA.is_infinite(fhA)) {
        std::cout << "No face located in A: " << from_exact(plane.to_3d(p)) << std::endl;
        //vout << " " << from_exact(plane.to_3d(p)) << std::endl;
      }
      if (fhA->info().id2 != std::make_pair(std::size_t(-1), std::size_t(-1))) {
        cit->info().idA2 = fhA->info().id2;
        // std::cerr << "A: " << fhA->info().id << std::endl;
        result.first += a;
      }
      else {
        std::cout << "Face in A is missing ID " << from_exact(plane.to_3d(p)) << std::endl;
        //vout << " " << from_exact(plane.to_3d(p)) << std::endl;
      }
      typename CDTplus::Face_handle fhB = cdtB.locate(p);
      if (cdtB.is_infinite(fhB)) {
        std::cout << "No face located in B: " << from_exact(plane.to_3d(p)) << std::endl;
        //vout << " " << from_exact(plane.to_3d(p)) << std::endl;
      }
      if (fhB->info().id2 != std::make_pair(std::size_t(-1), std::size_t(-1))) {
        cit->info().idB2 = fhB->info().id2;
        // std::cerr << "B: " << fhB->info().id << std::endl;
        result.second += a;
      }
      else {
        std::cout << "Face in B is missing ID " << from_exact(plane.to_3d(p)) << std::endl;
        //vout << " " << from_exact(plane.to_3d(p)) << std::endl;
      }
    }

    //vout.close();

    return result;
  }

  std::size_t check_cdt(CDTplus& cdt, const typename Intersection_kernel::Plane_3& plane) {
    std::size_t missing = 0;
    for (typename CDTplus::Finite_faces_iterator cit = cdt.finite_faces_begin(); cit != cdt.finite_faces_end(); ++cit) {
      if (cit->info().id2 == std::make_pair(std::size_t(-1), std::size_t(-1))) {
        std::cout << missing << ":";
        const std::string vfilename = std::to_string(missing) + "-missing-id.polylines.txt";
        std::ofstream vout(vfilename);
        vout.precision(20);
        vout << 4;
        for (std::size_t i = 0; i < 3; i++) {
          std::cout << " v(" << cit->vertex(i)->info().idx2.first << ", " << cit->vertex(i)->info().idx2.second << ")";
          vout << " " << plane.to_3d(cit->vertex(i)->point());
        }
        vout << " " << plane.to_3d(cit->vertex(0)->point()) << std::endl;
        std::cout << std::endl;
        vout << std::endl;
        vout.close();
        missing++;
      }
    }

    return missing;
  }

  std::size_t check_fusioned_cdt(CDTplus& cdt, const typename Intersection_kernel::Plane_3& plane) {
    std::size_t missing = 0;
    for (typename CDTplus::Finite_faces_iterator cit = cdt.finite_faces_begin(); cit != cdt.finite_faces_end(); ++cit) {
      if (cit->info().idA2 == std::make_pair(std::size_t(-1), std::size_t(-1)) || cit->info().idB2 == std::make_pair(std::size_t(-1), std::size_t(-1))) {
        std::cout << missing << ":";
        const std::string vfilename = std::to_string(missing) + "-missing-id.polylines.txt";
        std::ofstream vout(vfilename);
        vout.precision(20);
        vout << 4;
        for (std::size_t i = 0; i < 3; i++) {
          std::cout << " v(" << cit->vertex(i)->info().idx2.first << ", " << cit->vertex(i)->info().idx2.second << ")";
          vout << " " << plane.to_3d(cit->vertex(i)->point());
        }
        vout << " " << plane.to_3d(cit->vertex(0)->point()) << std::endl;
        std::cout << std::endl;
        vout << std::endl;
        vout.close();
        missing++;
      }
    }

    return missing;
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

  void collect_faces(std::size_t partition_idx, std::size_t sp_idx, std::vector<Index>& faces, typename Intersection_kernel::Plane_3& plane) {
    Sub_partition& p = m_partition_nodes[partition_idx];

    plane = p.m_data->support_plane(sp_idx).data().exact_plane;

    const std::vector<std::size_t>& f2sp = p.m_data->face_to_support_plane();

    for (std::size_t i = 0; i < f2sp.size(); i++)
      if (f2sp[i] == sp_idx)
        faces.push_back(std::make_pair(partition_idx, i));

/*
    for (std::size_t i = 0; i < p.m_data->volumes().size(); i++) {
      typename Data_structure::Volume_cell& v = p.m_data->volumes()[i];
      for (std::size_t j = 0; j < v.faces.size(); j++) {
        if (v.pfaces[j].first == sp_idx) {
          faces.push_back(std::make_pair(partition_idx, v.faces[j]));
        }
      }
    }*/
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
            std::cout << "neighbor is set partition: " << partition_idx << " volume: " << i << " face: " << j <<  " " << v.neighbors[j] << std::endl;
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

  void collect_faces(Octree_node node, std::size_t dimension, bool lower, std::vector<Index>& faces, typename Intersection_kernel::Plane_3 &plane) {
    // Collects boundary faces of node from its children.
    // dimension specifies the axis of the boundary face and lower determines if it is the lower of upper face of the cube on the axis.

    // Support plane indices:
    // xmin 4, xmax 2
    // ymin 1, ymax 3
    // zmin 0, zmax 5

    if (m_octree->is_leaf(node)) {
      // Mapping to partition is needed.
      std::size_t idx = m_node2partition[node];
      Sub_partition& partition = m_partition_nodes[m_node2partition[node]];
      From_exact from_exact;

      if (lower)
        switch (dimension) {
        case 0:
          collect_faces(idx, 4, faces, plane);
          break;
        case 1:
          collect_faces(idx, 1, faces, plane);
          break;
        case 2:
          collect_faces(idx, 0, faces, plane);
          break;
        }
      else
        switch (dimension) {
        case 0:
          collect_faces(idx, 2, faces, plane);
          break;
        case 1:
          collect_faces(idx, 3, faces, plane);
          break;
        case 2:
          collect_faces(idx, 5, faces, plane);
          break;
        }


      // Is Index as type for faces sufficient?
      // Partition/face id
      // Access to volumes via Data_structure->face_to_volumes()
      // However, the Data_structure::m_face2volumes needs to have std::pair<Index, Index> type to reference two volumes (pair<int, Index> would also be possible as only one volume can lie outside
      return;
    }
    else {
      typename Intersection_kernel::Plane_3 pl2, pl3, pl4;
      if (lower)
        switch (dimension) {
        case 0://0246
          collect_faces(m_octree->child(node, 0), dimension, true, faces, plane);
          collect_faces(m_octree->child(node, 2), dimension, true, faces, pl2);
          collect_faces(m_octree->child(node, 4), dimension, true, faces, pl3);
          collect_faces(m_octree->child(node, 6), dimension, true, faces, pl4);
          break;
        case 1://0145
          collect_faces(m_octree->child(node, 0), dimension, true, faces, plane);
          collect_faces(m_octree->child(node, 1), dimension, true, faces, pl2);
          collect_faces(m_octree->child(node, 4), dimension, true, faces, pl3);
          collect_faces(m_octree->child(node, 5), dimension, true, faces, pl4);
          break;
        case 2://0123
          collect_faces(m_octree->child(node, 0), dimension, true, faces, plane);
          collect_faces(m_octree->child(node, 1), dimension, true, faces, pl2);
          collect_faces(m_octree->child(node, 2), dimension, true, faces, pl3);
          collect_faces(m_octree->child(node, 3), dimension, true, faces, pl4);
          break;
        }
      else
        switch (dimension) {
        case 0://1357
          collect_faces(m_octree->child(node, 1), dimension, false, faces, plane);
          collect_faces(m_octree->child(node, 3), dimension, false, faces, pl2);
          collect_faces(m_octree->child(node, 5), dimension, false, faces, pl3);
          collect_faces(m_octree->child(node, 7), dimension, false, faces, pl4);
          break;
        case 1://3467
          collect_faces(m_octree->child(node, 2), dimension, false, faces, plane);
          collect_faces(m_octree->child(node, 3), dimension, false, faces, pl2);
          collect_faces(m_octree->child(node, 6), dimension, false, faces, pl3);
          collect_faces(m_octree->child(node, 7), dimension, false, faces, pl4);
          break;
        case 2://4567
          collect_faces(m_octree->child(node, 4), dimension, false, faces, plane);
          collect_faces(m_octree->child(node, 5), dimension, false, faces, pl2);
          collect_faces(m_octree->child(node, 6), dimension, false, faces, pl3);
          collect_faces(m_octree->child(node, 7), dimension, false, faces, pl4);
          break;
        }

      bool same = plane == pl2;
      same = (same && plane == pl3);
      same = (same && plane == pl4);
      if (!same) {
        std::cout << "collect_faces: different plane, node: " << node << " lower: " << lower << std::endl;
        From_exact from_exact;
        std::cout << from_exact(plane) << std::endl;
        std::cout << from_exact(pl2) << " child: " << m_octree->child(node, 4) << std::endl;
        std::cout << from_exact(pl3) << " child: " << m_octree->child(node, 6) << std::endl;
        std::cout << from_exact(pl4) << " child: " << m_octree->child(node, 7) << std::endl << std::endl;
      }
    }
  }

  void collect_opposing_faces(Octree_node node, std::size_t dimension, std::vector<Index>& lower, std::vector<Index>& upper, typename Intersection_kernel::Plane_3 &plane) {
    // Nothing to do for a leaf node.
    if (m_octree->is_leaf(node))
      return;

    typename Intersection_kernel::Plane_3 pl[7];
    switch (dimension) {
    case 0:
      collect_faces(m_octree->child(node, 0), dimension, false, lower, plane);
      collect_faces(m_octree->child(node, 2), dimension, false, lower, pl[0]);
      collect_faces(m_octree->child(node, 4), dimension, false, lower, pl[1]);
      collect_faces(m_octree->child(node, 6), dimension, false, lower, pl[2]);
      collect_faces(m_octree->child(node, 1), dimension, true, upper, pl[3]);
      collect_faces(m_octree->child(node, 3), dimension, true, upper, pl[4]);
      collect_faces(m_octree->child(node, 5), dimension, true, upper, pl[5]);
      collect_faces(m_octree->child(node, 7), dimension, true, upper, pl[6]);
      break;
    case 1:
      collect_faces(m_octree->child(node, 0), dimension, false, lower, plane);
      collect_faces(m_octree->child(node, 1), dimension, false, lower, pl[0]);
      collect_faces(m_octree->child(node, 4), dimension, false, lower, pl[1]);
      collect_faces(m_octree->child(node, 5), dimension, false, lower, pl[2]);
      collect_faces(m_octree->child(node, 3), dimension, true, upper, pl[3]);
      collect_faces(m_octree->child(node, 2), dimension, true, upper, pl[4]);
      collect_faces(m_octree->child(node, 6), dimension, true, upper, pl[5]);
      collect_faces(m_octree->child(node, 7), dimension, true, upper, pl[6]);
      break;
    case 2:
      collect_faces(m_octree->child(node, 0), dimension, false, lower, plane);
      collect_faces(m_octree->child(node, 1), dimension, false, lower, pl[0]);
      collect_faces(m_octree->child(node, 2), dimension, false, lower, pl[1]);
      collect_faces(m_octree->child(node, 3), dimension, false, lower, pl[2]);
      collect_faces(m_octree->child(node, 4), dimension, true, upper, pl[3]);
      collect_faces(m_octree->child(node, 5), dimension, true, upper, pl[4]);
      collect_faces(m_octree->child(node, 6), dimension, true, upper, pl[5]);
      collect_faces(m_octree->child(node, 7), dimension, true, upper, pl[6]);
      break;
    }

    From_exact from_exact;
    //std::cout << from_exact(plane) << std::endl;

    bool same = true;
    for (std::size_t i = 0; i < 3; i++)
      same = (same && plane == pl[i]);

    for (std::size_t i = 3; i < 7; i++)
      same = (same && plane.opposite() == pl[i]);

    if (!same) {
      std::cout << "collect_opposing_faces: different plane, node: " << node << std::endl;
      std::cout << from_exact(plane) << std::endl;
      for (std::size_t i = 0; i < 3; i++)
        std::cout << from_exact(pl[i]) << std::endl;
      for (std::size_t i = 3; i < 7; i++)
        std::cout << from_exact(pl[i].opposite()) << std::endl;
      bool diff = (plane.b() == pl[6].opposite().b());
      std::cout << diff << std::endl;
      std::cout << std::endl;
    }
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

      CDTplus lowerCDT, upperCDT;
      double lower_area = build_cdt(lowerCDT, vertices_lower, faces_lower, plane);
      double upper_area = build_cdt(upperCDT, vertices_upper, faces_upper, plane);

      // Collect Plane_3 from support planes (two vectors, one for each other)
      std::vector<std::size_t> planes_lower, planes_upper;
      collect_planes(lower_idx, lower_sp, planes_lower);
      collect_planes(upper_idx, upper_sp, planes_upper);

      // Remove common planes
      auto lower_it = planes_lower.begin();
      auto upper_it = planes_upper.begin();
      while (lower_it != planes_lower.end() && upper_it != planes_upper.end()) {
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

      //CDTplus lowerCDT, upperCDT;
      lower_area = build_cdt(lowerCDT, vertices_lower, faces_lower, plane);
      upper_area = build_cdt(upperCDT, vertices_upper, faces_upper, plane);

      CDTplus combined;
      std::pair<double, double> areas = overlay(combined, lowerCDT, upperCDT, plane);

      for (std::size_t i = 0; i < faces_lower.size(); i++) {
        typename Data_structure::Volume_cell& v = lower.m_data->volumes()[face_idx_lower[i].first];

        std::vector<typename Intersection_kernel::Point_3> pts;
        pts.reserve(faces_lower[i].size());
        for (std::size_t idx : faces_lower[i])
          pts.push_back(vertices_lower[idx]);

        typename Intersection_kernel::Point_3 c = CGAL::centroid(pts.begin(), pts.end(), CGAL::Dimension_tag<0>());
        Point_3 c_inexact = from_exact(c);

        typename CDTplus::Face_handle neighbor = upperCDT.locate(plane.to_2d(c));
        if (neighbor->info().id < faces_upper.size()) {
          //std::cout << "index " << i << " of lower set to " << face_idx_upper[neighbor->info().id].first << std::endl;
          v.neighbors[face_idx_lower[i].second] = face_idx_upper[neighbor->info().id].first;
        }
        else std::cout << "neighbor of face " << i << " of lower has neighbor " << face_idx_upper[neighbor->info().id].first << " in upper" << std::endl;
      }

      //check_faces(lower_idx, lower_sp);

      for (std::size_t i = 0; i < faces_upper.size(); i++) {
        typename Data_structure::Volume_cell& v = upper.m_data->volumes()[face_idx_upper[i].first];

        std::vector<typename Intersection_kernel::Point_3> pts;
        pts.reserve(faces_upper[i].size());
        for (std::size_t idx : faces_upper[i])
          pts.push_back(vertices_upper[idx]);

        typename Intersection_kernel::Point_3 c = CGAL::centroid(pts.begin(), pts.end(), CGAL::Dimension_tag<0>());

        typename CDTplus::Face_handle neighbor = lowerCDT.locate(plane.to_2d(c));
        if (neighbor->info().id < faces_lower.size()) {
          //std::cout << "index " << i << " of upper set to " << face_idx_lower[neighbor->info().id].first << std::endl;
          v.neighbors[face_idx_upper[i].second] = face_idx_lower[neighbor->info().id].first;
        }
        else std::cout << "neighbor of face " << i << " of upper has neighbor " << face_idx_lower[neighbor->info().id].first << " in upper" << std::endl;
      }

      //check_faces(upper_idx, upper_sp);
    }
  }

  bool same_face(const Face_handle& a, const Face_handle& b) const {
    return (b->info().idA2 == a->info().idA2 && b->info().idB2 == a->info().idB2);
  }

  void dump_face(const Face_handle& f, const std::string& filename) {
    From_exact from_exact;
    std::ofstream vout(filename);
    vout.precision(20);
    vout << "4 ";
    vout << " " << from_exact(f->vertex(0)->info().point_3);
    vout << " " << from_exact(f->vertex(1)->info().point_3);
    vout << " " << from_exact(f->vertex(2)->info().point_3);
    vout << " " << from_exact(f->vertex(0)->info().point_3);
    vout << std::endl;
    vout.close();
  }

  void dump_point(const Vertex_handle& v, const std::string& filename) {
    From_exact from_exact;
    std::ofstream vout3(filename);
    vout3.precision(20);
    vout3 << " " << from_exact(v->info().point_3);
    vout3 << std::endl;
    vout3.close();
  }

  void set_face(const Index& f, const Index& other, std::set<Index>& replaced, const std::vector<Vertex_handle>& polygon) {
    From_exact from_exact;
    auto pair = replaced.insert(f);
    std::size_t idx;
    assert(m_partition_nodes[f.first].face_neighbors[f.second].first.first == f.first);
    std::size_t vol_idx = m_partition_nodes[f.first].face_neighbors[f.second].first.second;
    if (!pair.second) {
      // New face has a new index
      idx = m_partition_nodes[f.first].face2vertices.size();
      // Add face into vector
      m_partition_nodes[f.first].face2vertices.push_back(std::vector<Index>());
      m_partition_nodes[f.first].m_data->face_is_part_of_input_polygon().push_back(m_partition_nodes[f.first].m_data->face_is_part_of_input_polygon()[f.second]);
      // Add face index into volume
      m_partition_nodes[f.first].m_data->volumes()[vol_idx].faces.push_back(idx);
      // Copy neighbor from already existing face
      m_partition_nodes[f.first].face_neighbors.push_back(m_partition_nodes[f.first].face_neighbors[f.second]);
      m_partition_nodes[f.first].m_data->face_to_support_plane().push_back(m_partition_nodes[f.first].m_data->face_to_support_plane()[f.second]);
    }
    else {
      idx = f.second;
      // All boundary faces should have a negative second neighbor.
      assert(m_partition_nodes[f.first].face_neighbors[idx].second.second >= std::size_t(-6));
    }
    std::vector<Index>& vertices = m_partition_nodes[f.first].face2vertices[idx];
    // First neighbor of other face should point to the inside volume in the other partition and thus cannot be negative
    assert(m_partition_nodes[other.first].face_neighbors[other.second].first.second < std::size_t(-6));
    m_partition_nodes[f.first].face_neighbors[idx].second = m_partition_nodes[other.first].face_neighbors[other.second].first;
    vertices.resize(polygon.size());
    for (std::size_t i = 0; i < polygon.size(); i++) {
      VI& vi = polygon[i]->info();
      // Is this check actually meaningless as partition indices now start at 0?
      // Check whether they are initialized as 0 and where it is used as indicator for something.
/*
      if (vi.idA2.first == 0 || vi.idB2.first == 0) {
        std::cout << "invalid vertex id" << std::endl;
      }*/

      if (vi.idA2.first != std::size_t(-1))
        vertices[i] = vi.idA2;
      else if (vi.idB2.first != std::size_t(-1))
        vertices[i] = vi.idB2;
      else {
        std::size_t vidx = m_partition_nodes[f.first].m_data->vertices().size();
        m_partition_nodes[f.first].m_data->vertices().push_back(from_exact(vi.point_3));
        m_partition_nodes[f.first].m_data->exact_vertices().push_back(vi.point_3);
        vertices[i] = vi.idA2 = std::make_pair(f.first, vidx);
        // Prevent T-junctions here!
        // Check adjacent volumes, adjacent partitions
        // There is no edge connectivity information
      }
    }
  }

  void adapt_faces(const CDTplus& cdt, std::vector<Index>& a, std::vector<Index>& b, typename Intersection_kernel::Plane_3& plane) {
    std::set<Index> replacedA, replacedB;
    From_exact from_exact;

    std::size_t extracted = 0;
    for (typename CDTplus::Face_handle fh : cdt.finite_face_handles()) {
      // when extracting each face, I only need to insert vertices, that don't exist on either side. Otherwise, I can just reference the vertex in the other partition.
      // using cit->info().id2.first as visited flag. -1 means not visited
      if (fh->info().id2.first != -1)
        continue;

      // 4 different cases: no border edge, 1, 2 or 3
      // Check 1, 2, 3
      // if first is not, continue
      // if first move in other direction? search start

      // Easier approach, don't make a list of edges, but a list of vertices
      // Find first pair of vertices, then just loop around last vertex using Face_circulator
      // -> Triangulation only has vertices and faces, no easy way to loop over edges

      std::vector<Vertex_handle> face;

      for (std::size_t i = 0; i < 3; i++)
        if (cdt.is_infinite(fh->neighbor(i)) || !same_face(fh, fh->neighbor(i))) {
          face.push_back(fh->vertex((i + 2) % 3));
          face.push_back(fh->vertex((i + 1) % 3));
          break;
        }

      // No border edge?
      if (face.empty())
        continue;
      else {
        //dump_point(face.back(), "last.xyz");
        Face_handle last = fh;

        // Mark seed face as segmented
        fh->info().id2.first = extracted;

        // edge is pair<Face_handle, int (vertex)>
        while (face.front() != face.back()) {
          auto eit = cdt.incident_edges(face.back(), last);
          // Skip the first edge as it always starts with the edge itself.
          //eit++;
/*
          auto eit2 = eit;
          for (std::size_t i = 0; i < 10; i++) {
            dump_point(eit2->first->vertex(eit2->second), std::to_string(i) + "p.xyz");
            dump_face(eit2->first, std::to_string(i) + "tri.polylines.txt");
            std::cout << i << " same: " << same_face(last, eit2->first->neighbor((eit2->second + 1) % 3)) << std::endl;
            eit2++;
          }*/
          auto first = eit;
          Point_3 p = from_exact(eit->first->vertex(eit->second)->info().point_3);

          assert(!cdt.is_infinite(eit->first));
          do {
            // Export tri

            //dump_point(eit->first->vertex(eit->second), "p.xyz");
            //dump_face(eit->first, "tri.polylines.txt");

            // Is the current edge to the infinite vertex?
            if (cdt.is_infinite(eit->first->neighbor((eit->second + 1) % 3))) {
              eit++;
              continue;
            }

            bool infinite = cdt.is_infinite(eit->first);

            /*            if (!infinite)
                          dump_face(eit->first, "neighbor.polylines.txt");*/

            if (infinite || !same_face(last, eit->first)) {
              last = eit->first->neighbor((eit->second + 1) % 3);
              last->info().id2.first = extracted;
              face.push_back(eit->first->vertex(eit->second));

              break;
            }
            eit++;
            assert(eit != first);
          } while (eit != first);
          // If last vertex is equal to first vertex, stop
          // Take last vertex and face
          // First find index of vertex in that face
          // Check if opposite face of next edge, if not same, add next vertex and reloop
          // if not, check next face

          assert(face.size() < 100);
        }

        // The last vertex is equal to the first one, so it should be removed.
        face.pop_back();

        // face ids in partitions are stored in fh->info
        ID& id = fh->info();
        set_face(id.idA2, id.idB2, replacedA, face);
        set_face(id.idB2, id.idA2, replacedB, face);
      }

      // Checking for border edges. If opposite faces do not exist or don't have the same indices, the edge belongs to a new face.
      // cit->neighbor(i) is the face opposite of vertex(i), meaning on the other side of the edge between vertex((i+1)%3) and vertex((i+2)%3)
    }
  }

  void make_conformal(std::vector<Index>& a, std::vector<Index>& b, typename Intersection_kernel::Plane_3 &plane) {
    // partition ids are in a[0].first and b[0].first
    // volume and face in volume ids are not available
    // there is face2volume and one of those volume indices will be an outside volume, e.g. std::size_t(-1) to std::size_t(-6)

    // Indices in a and b come from different partitions. Each face only has vertices from the same partition
    // Constraints in the cdt should have matching vertices and edges from different partitions -> opportunity to match vertices and faces between partitions

    // buildCDT needs only Index and exact_vertices for the points and Index for faces

    std::unordered_map<std::size_t, std::vector<Index> > a_sets, b_sets;
    for (const Index& i : a) {
      a_sets[i.first].push_back(i);
      if (m_partition_nodes[i.first].m_data->face_is_part_of_input_polygon()[i.second])
        std::cout << "(" << i.first << ", " << i.second << ") is part of input polygon" << std::endl;
    }
    for (const Index& i : b) {
      b_sets[i.first].push_back(i);
      if (m_partition_nodes[i.first].m_data->face_is_part_of_input_polygon()[i.second])
        std::cout << "(" << i.first << ", " << i.second << ") is part of input polygon" << std::endl;
    }

    std::vector<CDTplus> a_cdts(a_sets.size()), b_cdts(b_sets.size());

    Index g(-1, -1);
    std::size_t newpts = 0;
    From_exact from_exact;
    Plane_3 pl = from_exact(plane);

    std::size_t idx = 0;
    for (auto& p : a_sets) {
      build_cdt(a_cdts[idx], p.second, plane);
      newpts = 0;
      for (Vertex_handle v : a_cdts[idx].finite_vertex_handles()) {
        if (v->info().idA2 == g)
          newpts++;
      }

      if (newpts > 0)
        std::cout << newpts << " vertices without references found in a_cdts" << idx << std::endl;

      if (check_cdt(a_cdts[idx], plane) != 0)
        std::cout << "lower " << p.first << ": " << p.second.size() << " " << a_cdts[idx].number_of_faces() << " with " << check_cdt(a_cdts[idx], plane) << " missing ids" << std::endl;
      idx++;
    }

    idx = 0;
    for (auto& p : b_sets) {
      build_cdt(b_cdts[idx], p.second, plane);

      newpts = 0;
      for (Vertex_handle v : b_cdts[idx].finite_vertex_handles()) {
        if (v->info().idA2 == g)
          newpts++;
      }

      if (newpts > 0)
        std::cout << newpts << " vertices without references found in b_cdts" << idx << std::endl;

      if (check_cdt(b_cdts[idx], plane) != 0)
        std::cout << "upper " << p.first << ": " << p.second.size() << " " << b_cdts[idx].number_of_faces() << " with " << check_cdt(b_cdts[idx], plane) << " missing ids" << std::endl;
      idx++;
    }

    CDTplus cdtA, cdtB, cdtC;
    build_cdt(cdtA, a_cdts, plane);
    std::size_t missing = check_cdt(cdtA, plane);
    if (missing > 0)
      std::cout << "lower: " << a.size() << " " << cdtA.number_of_faces() << " faces " << cdtA.number_of_vertices() << " vertices with " << missing << " missing ids" << std::endl;

/*
    std::ofstream vout("cdtA.polylines.txt");
    vout.precision(20);
    for (typename CDTplus::Face_handle fh : cdtA.finite_face_handles()) {
      vout << "4 ";
      vout << " " << from_exact(fh->vertex(0)->info().point_3);
      vout << " " << from_exact(fh->vertex(1)->info().point_3);
      vout << " " << from_exact(fh->vertex(2)->info().point_3);
      vout << " " << from_exact(fh->vertex(0)->info().point_3);
      vout << std::endl;
    }
    vout << std::endl;
    vout.close();*/

/*
    for (Vertex_handle v : cdtA.finite_vertex_handles()) {
      if (v->info().idA2 == g && v->info().idB2 == g)
        newpts++;
    }

    std::cout << newpts << " vertices without references found in cdtA" << std::endl;*/

    build_cdt(cdtB, b_cdts, plane);
    missing = check_cdt(cdtB, plane);
    if (missing > 0)
      std::cout << "upper: " << b.size() << " " << cdtB.number_of_faces() << " faces " << cdtB.number_of_vertices() << " vertices with " << missing << " missing ids" << std::endl;

/*
    std::ofstream vout2("cdtB.polylines.txt");
    vout2.precision(20);
    for (typename CDTplus::Face_handle fh : cdtB.finite_face_handles()) {
      vout2 << "4 ";
      vout2 << " " << from_exact(fh->vertex(0)->info().point_3);
      vout2 << " " << from_exact(fh->vertex(1)->info().point_3);
      vout2 << " " << from_exact(fh->vertex(2)->info().point_3);
      vout2 << " " << from_exact(fh->vertex(0)->info().point_3);
      vout2 << std::endl;
    }
    vout2 << std::endl;
    vout2.close();*/

/*
    newpts = 0;
    for (Vertex_handle v : cdtB.finite_vertex_handles()) {
      if (v->info().idA2 == g && v->info().idB2 == g)
        newpts++;
    }

    std::cout << newpts << " vertices without references found in cdtB" << std::endl;*/

    overlay(cdtC, cdtA, cdtB, plane);
    //std::cout << "overlay: " << cdtC.number_of_faces() << " faces " << cdtC.number_of_vertices() << " vertices" << std::endl;

    adapt_faces(cdtC, a, b, plane);

    // Is there linkage between the cdts? I could create a map of vertex Index to cdt vertices
    // I can create an unordered map from face Index to vector of cdt_face iterator

    // Each input face can be split into several faces
    // Updating the neighbor volumes does not seem difficult but not completely trivial either as it has to be done after the face extraction (due to the neighbors array in volumes)
    // -> each face extracted has the same volume ids (or the same face ids on both sides)

    // Walk around the edges to identify faces
    // - segment by identical face ids on both sides
    // How to keep track of the face vector in volumes? I can change the first face in place

    // How to identify edges that are split? Can I add a property on edges to mark that they have been added? Seems difficult, because the vertex can be part of plenty new edges.
    // Can it? The overlay is a fusion of 2 cdts, so if there are more than two edges intersecting in a vertex, there were already two edges intersecting in one of the cdts
    // So no, each new vertex can only be part of 2 edges

    // Adjusting edges part of faces that are not part of the splitting plane is basically a function of split_edge(Index_head, Index_tail, new_mid_vertex)
    // Identifying the faces based on the vertices seems costly. PEdge to PFaces exists, check for PFace to volume
    // // -> does not work for newly inserted edges! Newly inserted edges do not have PEdge!
    // Otherwise it is possible to find adjacent volumes based on the participating faces. However, an edge on the boundary can be part of many faces/volumes

    // Approach:
    // Loop over finite faces of fusioned cdt
    // check if face was already handled (-> reuse field in face info?)
    // check if face is on border to face of another index pair
    //  start face extraction
    //   follow border and build up polygon vector
    //   check if there is a vertex index in the vertex info, if not insert vertex into partition.data_structure and update
    //   create map for boundary vertices correspondences
    //  check if replace face in data structure or create new one (set which contains replaced ones?)
  }

  void make_conformal(Octree_node node) {
    // Nothing to do for a leaf node.
    if (m_octree->is_leaf(node))
      return;

    // Make childs conformal
    for (std::size_t i = 0; i < 8; i++) {
      make_conformal(m_octree->child(node, i));
    }

    // Make itself conformal
    // Get faces between child nodes
    // do in inverse dimension order (like inverse splitting order, start by 2 or 1 and walk down to 0)
    // follow cdt approach in split_octree

    // Order of children?
    // x, y, z planes can be merged independently
    for (std::size_t dim = 0; dim < 3; dim++) {
      std::vector<Index> lower, upper;
      typename Intersection_kernel::Plane_3 plane;

      collect_opposing_faces(node, dim, lower, upper, plane);

      make_conformal(lower, upper, plane);

      lower.clear();
      upper.clear();
      collect_opposing_faces(node, dim, lower, upper, plane);

      for (std::size_t i = 0; i < lower.size(); i++) {
        auto n = neighbors(lower[i]);
        assert(n.first >= 0 && n.second >= 0);
      }

      for (std::size_t i = 0; i < upper.size(); i++) {
        auto n = neighbors(upper[i]);
        assert(n.first >= 0 && n.second >= 0);
      }
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

  void split_octree() {
    // Octree creation for sub partition
    std::size_t count = 0;
    for (const auto& p : m_input_polygons)
      count += p.size();

    m_points.clear();
    m_points.reserve(count);
    m_polygons.reserve(m_input_polygons.size());

    for (const auto& p : m_input_polygons) {
      std::size_t idx = m_points.size();
      std::copy(p.begin(), p.end(), std::back_inserter(m_points));
      m_polygons.push_back(std::vector<std::size_t>(p.size()));
      std::iota(m_polygons.back().begin(), m_polygons.back().end(), idx);
    }

    m_octree = std::make_unique<Octree>(CGAL::Orthtree_traits_polygons<Kernel>(m_points, m_polygons, m_parameters.bbox_dilation_ratio));
    m_octree->refine(0, 40);

    /*
    // Collect all the leaf nodes
    std::queue<Node_index> leaf_nodes;
    for (Node_index leaf: traverse(Orthtrees::Leaves_traversal<Self>(*this))) {
      leaf_nodes.push(leaf);
    }
    */

    std::size_t leaf_count = 0;
    std::size_t max_count = 0;

    for (Octree::Node_index node : m_octree->traverse(CGAL::Orthtrees::Leaves_traversal<Octree>(*m_octree))) {
      if (m_octree->is_leaf(node))
        leaf_count++;
      else
        std::cout << "Leaves_traversal traverses non-leaves" << std::endl;
      max_count = (std::max<std::size_t>)(max_count, node);
    }

    m_partition_nodes.resize(leaf_count);

    m_node2partition.resize(max_count + 1, std::size_t(-1));

    std::size_t idx = 0;
    for (Octree::Node_index node : m_octree->traverse(CGAL::Orthtrees::Leaves_traversal<Octree>(*m_octree)))
      if (m_octree->is_leaf(node)) {
        // Creating bounding box
        CGAL::Iso_cuboid_3<Kernel> box = m_octree->bbox(node);
        m_partition_nodes[idx].bbox[0] = typename Intersection_kernel::Point_3(box.xmin(), box.ymin(), box.zmin());
        m_partition_nodes[idx].bbox[1] = typename Intersection_kernel::Point_3(box.xmax(), box.ymin(), box.zmin());
        m_partition_nodes[idx].bbox[2] = typename Intersection_kernel::Point_3(box.xmax(), box.ymax(), box.zmin());
        m_partition_nodes[idx].bbox[3] = typename Intersection_kernel::Point_3(box.xmin(), box.ymax(), box.zmin());
        m_partition_nodes[idx].bbox[4] = typename Intersection_kernel::Point_3(box.xmin(), box.ymax(), box.zmax());
        m_partition_nodes[idx].bbox[5] = typename Intersection_kernel::Point_3(box.xmin(), box.ymin(), box.zmax());
        m_partition_nodes[idx].bbox[6] = typename Intersection_kernel::Point_3(box.xmax(), box.ymin(), box.zmax());
        m_partition_nodes[idx].bbox[7] = typename Intersection_kernel::Point_3(box.xmax(), box.ymax(), box.zmax());

/*
        auto bbox = m_octree->bbox(i);
        m_partition_nodes[idx].bbox[0] = Point_3(bbox.xmin(), bbox.ymin(), bbox.zmin());
        m_partition_nodes[idx].bbox[1] = Point_3(bbox.xmax(), bbox.ymin(), bbox.zmin());
        m_partition_nodes[idx].bbox[2] = Point_3(bbox.xmax(), bbox.ymax(), bbox.zmin());
        m_partition_nodes[idx].bbox[3] = Point_3(bbox.xmin(), bbox.ymax(), bbox.zmin());
        m_partition_nodes[idx].bbox[4] = Point_3(bbox.xmin(), bbox.ymax(), bbox.zmax());
        m_partition_nodes[idx].bbox[5] = Point_3(bbox.xmin(), bbox.ymin(), bbox.zmax());
        m_partition_nodes[idx].bbox[6] = Point_3(bbox.xmax(), bbox.ymin(), bbox.zmax());
        m_partition_nodes[idx].bbox[7] = Point_3(bbox.xmax(), bbox.ymax(), bbox.zmax());*/

        // Get consistent Plane_3 from Octree to generate exact planes

        auto polys = m_octree->data(node);
        for (std::size_t j = 0; j < polys.size(); j++) {
          m_partition_nodes[idx].input_polygons.push_back(polys[j].first);
          m_partition_nodes[idx].m_input_planes.push_back(m_input_planes[polys[j].first]);
        }

        m_partition_nodes[idx].clipped_polygons.resize(polys.size());
        for (std::size_t i = 0; i < polys.size(); i++) {
          m_partition_nodes[idx].clipped_polygons[i].resize(polys[i].second.size());
          for (std::size_t j = 0; j < polys[i].second.size(); j++)
            m_partition_nodes[idx].clipped_polygons[i][j] = polys[i].second[j];
        }

        // set node index
        m_partition_nodes[idx].node = node;
        m_node2partition[node] = idx;

        if (m_parameters.debug) {
          const std::string vfilename = std::to_string(idx) + "-box.polylines.txt";
          std::ofstream vout(vfilename);
          vout.precision(20);
          // zmin side
          vout << 5;
          vout << " " << m_partition_nodes[idx].bbox[0];
          vout << " " << m_partition_nodes[idx].bbox[1];
          vout << " " << m_partition_nodes[idx].bbox[2];
          vout << " " << m_partition_nodes[idx].bbox[3];
          vout << " " << m_partition_nodes[idx].bbox[0];
          // zmax side
          vout << std::endl << 5;
          vout << " " << m_partition_nodes[idx].bbox[4];
          vout << " " << m_partition_nodes[idx].bbox[5];
          vout << " " << m_partition_nodes[idx].bbox[6];
          vout << " " << m_partition_nodes[idx].bbox[7];
          vout << " " << m_partition_nodes[idx].bbox[4];
          // 4 edges between zmin and zmax
          vout << std::endl << 2;
          vout << " " << m_partition_nodes[idx].bbox[0];
          vout << " " << m_partition_nodes[idx].bbox[5];
          vout << std::endl << 2;
          vout << " " << m_partition_nodes[idx].bbox[1];
          vout << " " << m_partition_nodes[idx].bbox[6];
          vout << std::endl << 2;
          vout << " " << m_partition_nodes[idx].bbox[2];
          vout << " " << m_partition_nodes[idx].bbox[7];
          vout << std::endl << 2;
          vout << " " << m_partition_nodes[idx].bbox[3];
          vout << " " << m_partition_nodes[idx].bbox[4];
          vout << std::endl;
          vout.close();

          KSR_3::dump_polygons(m_partition_nodes[idx].clipped_polygons, std::to_string(idx) + "-polys.ply");
        }
        idx++;
      }

    std::cout << "input split into " << m_partition_nodes.size() << " partitions" << std::endl;
  }

  bool within_tolerance(const Plane_3& p1, const Point_2 &c1, const Plane_3& p2, const Point_2& c2) const {
    using FT = typename GeomTraits::FT;

    const auto va = p1.orthogonal_vector();
    const auto vb = p2.orthogonal_vector();

    // Are the planes parallel?
    // const FT vtol = KSR::vector_tolerance<FT>();
    // const FT aval = CGAL::abs(va * vb);

    // std::cout << "aval: " << aval << " : " << vtol << std::endl;
    // if (aval < vtol) {
    //   return false;
    // }

    FT aval = approximate_angle(va, vb);
    CGAL_assertion(aval >= FT(0) && aval <= FT(180));
    if (aval >= FT(90))
      aval = FT(180) - aval;

    if (aval >= m_parameters.angle_tolerance) {
      return false;
    }

    const auto pa1 = p1.to_3d(c1);
    const auto pb1 = p2.projection(pa1);
    const auto pb2 = p2.to_3d(c2);
    const auto pa2 = p1.projection(pb2);

    const FT bval1 = KSR::distance(pa1, pb1);
    const FT bval2 = KSR::distance(pa2, pb2);
    const FT bval = (CGAL::max)(bval1, bval2);
    CGAL_assertion(bval >= FT(0));

    if (bval >= m_parameters.distance_tolerance)
      return false;

    return true;
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
