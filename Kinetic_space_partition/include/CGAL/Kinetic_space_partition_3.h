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

#ifndef CGAL_KINETIC_SPACE_PARTITION_3_H
#define CGAL_KINETIC_SPACE_PARTITION_3_H

#include <CGAL/license/Kinetic_space_partition.h>

// Boost includes.
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Named_function_parameters.h>
//#include <boost/filesystem.hpp>

#include <algorithm>
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
#include <CGAL/KSP/utils.h>
#include <CGAL/KSP/debug.h>
#include <CGAL/KSP/parameters.h>

#include <CGAL/KSP_3/Data_structure.h>
#include <CGAL/KSP_3/Initializer.h>
#include <CGAL/KSP_3/FacePropagation.h>
#include <CGAL/KSP_3/Finalizer.h>

#include <CGAL/Octree.h>
#include <CGAL/Orthtree_traits_polygons.h>

namespace CGAL {

/*!
* \ingroup PkgKineticSpacePartitionRef
  \brief creates the kinetic partition of the bounding box of the polygons given as input data. The kinetic partition can either be initialized
  by using the default constructor \link CGAL::Kinetic_space_partition_3::Kinetic_space_partition_3() `Kinetic_space_partition_3()`\endlink, `insert()` to provide input data and `initialize()` to prepare the partition or by using the constructor with input parameters.

  \tparam GeomTraits
    must be a model of `KineticSpacePartitionTraits_3`.

  \tparam IntersectionTraits
    must be a model of `Kernel` using exact computations. Defaults to `CGAL::Exact_predicates_exact_constructions_kernel`.
*/
template<typename GeomTraits, typename IntersectionTraits = CGAL::Exact_predicates_exact_constructions_kernel>
class Kinetic_space_partition_3 {

public:
  using Kernel = GeomTraits;
  using Intersection_kernel = IntersectionTraits;

  using Point_3 = typename Kernel::Point_3;

  using Index = std::pair<std::size_t, std::size_t>;

   /*!
   identifies the support of a face in the exported linear cell complex, which is either an input polygon, identified by a non-negative number, a side of the bounding box in the rotated coordinate system or a face of the octree used to partition the scene.
   */
  enum Face_support : int {
    ZMIN        = -1,
    YMIN        = -2,
    XMAX        = -3,
    YMAX        = -4,
    XMIN        = -5,
    ZMAX        = -6,
    OCTREE_FACE = -7,
  };

  /*!
  \brief this class provides a minimal model of `KineticLCCItems`. It adds attributes to faces and volumes and defines the use of index-based `LinearCellComplex`.
  */
  class Linear_cell_complex_min_items {
  public:
    typedef CGAL::Tag_true Use_index;
    typedef std::uint32_t Index_type;

    struct Face_attribute {
      Face_support input_polygon_index; // Non-negative numbers represent the index of the input polygon. Negative numbers correspond to the values defined in the enum `Face_support`.
      typename Intersection_kernel::Plane_3 plane;
      bool part_of_initial_polygon;
    };

    struct Volume_attribute {
      typename Intersection_kernel::Point_3 barycenter;
      std::size_t volume_id;
    };

    template<class LCC>
    struct Dart_wrapper {
      typedef CGAL::Cell_attribute_with_point< LCC, void > Vertex_cell_attribute;
      typedef CGAL::Cell_attribute< LCC, Face_attribute > Face_cell_attribute;
      typedef CGAL::Cell_attribute< LCC, Volume_attribute > Volume_cell_attribute;

      typedef std::tuple<Vertex_cell_attribute, void, Face_cell_attribute, Volume_cell_attribute> Attributes;
    };
  };

private:
  using FT = typename Kernel::FT;
  using Point_2 = typename Kernel::Point_2;
  using Vector_2 = typename Kernel::Vector_2;
  using Vector_3 = typename Kernel::Vector_3;
  using Plane_3 = typename Kernel::Plane_3;
  using Line_3 = typename Kernel::Line_3;
  using Line_2 = typename Kernel::Line_2;
  using Triangle_2 = typename Kernel::Triangle_2;
  using Transform_3 = CGAL::Aff_transformation_3<Kernel>;

  using Data_structure = KSP_3::internal::Data_structure<Kernel, Intersection_kernel>;

  using IVertex = typename Data_structure::IVertex;
  using IEdge   = typename Data_structure::IEdge;

  using From_exact = typename CGAL::Cartesian_converter<Intersection_kernel, Kernel>;
  using To_exact = typename CGAL::Cartesian_converter<Kernel, Intersection_kernel>;

  using Initializer = KSP_3::internal::Initializer<Kernel, Intersection_kernel>;
  using Propagation = KSP_3::internal::FacePropagation<Kernel, Intersection_kernel>;
  using Finalizer   = KSP_3::internal::Finalizer<Kernel, Intersection_kernel>;

  using Polygon_mesh = CGAL::Surface_mesh<Point_3>;
  using Timer        = CGAL::Real_timer;
  using Parameters = KSP::internal::Parameters_3<FT>;

  using Octree = CGAL::Orthtree<CGAL::Orthtree_traits_polygons<Kernel> >;
  using Octree_node = typename Octree::Node_index;

  struct VI
  {
    VI()
      : idA2(-1, -1), idB2(-1, -1), input(false)
    {}

    void set_point(const typename Intersection_kernel::Point_3& p) {
      point_3 = p;
      input = true;
    }

    typename Intersection_kernel::Point_3 point_3;
    std::set<Index> adjacent;
    Index idA2, idB2;
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
    int id, idA, idB;
    Index id2, idA2, idB2;
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
  CGAL::Aff_transformation_3<Intersection_kernel> m_transform;
  std::vector<Sub_partition> m_partition_nodes; // Tree of partitions.
  std::vector<std::size_t> m_partitions; // Contains the indices of the leaf nodes, the actual partitions to be calculated.
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

  std::set<Index> duplicates;

public:
  /// \name Initialization
  /// @{
  /*!
  \brief constructs an empty kinetic space partition object. Use `insert()` afterwards to insert polygons into the partition and `initialize()` to initialize the partition.

  \tparam NamedParameters
  a sequence of \ref bgl_namedparameters "Named Parameters"

  \param np
  a sequence of \ref bgl_namedparameters "Named Parameters"
  among the ones listed below

  \cgalNamedParamsBegin
    \cgalParamNBegin{verbose}
      \cgalParamDescription{Write timing and internal information to `std::cout`.}
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
  Kinetic_space_partition_3(
    const NamedParameters& np = CGAL::parameters::default_values()) :
    m_parameters(
      parameters::choose_parameter(parameters::get_parameter(np, internal_np::verbose), false),
      parameters::choose_parameter(parameters::get_parameter(np, internal_np::debug), false)), // use true here to export all steps
    m_input2regularized() {}

  /*!
  \brief constructs a kinetic space partition object and initializes it.

  \tparam PointRange
  must be a model of `ConstRange` whose iterator type is `RandomAccessIterator` and whose value type is Point_3.

  \tparam PolygonRange
  contains index ranges to form polygons by providing indices into PointRange

  \tparam NamedParameters
  a sequence of \ref bgl_namedparameters "Named Parameters"

  \param points
  an instance of `PointRange` with 3D points and corresponding 3D normal vectors

  \param polygons
  a range of non-coplanar polygons defined by a range of indices into `points`

  \param np
  a sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

  \cgalNamedParamsBegin
     \cgalParamNBegin{point_map}
       \cgalParamDescription{a property map associating points to the elements of the input range PointRange `points`.}
       \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type of the iterator of `PointRange` and whose value type is `GeomTraits::Point_3`}
       \cgalParamDefault{`CGAL::Identity_property_map<GeomTraits::Point_3>`}
     \cgalParamNEnd
    \cgalParamNBegin{debug}
      \cgalParamDescription{Export of intermediate results.}
      \cgalParamType{bool}
      \cgalParamDefault{false}
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
    \cgalParamNBegin{max_octree_depth}
      \cgalParamDescription{The maximal depth of the octree for subdividing the kinetic partition before initialization.}
      \cgalParamType{std::size_t}
      \cgalParamDefault{3}
    \cgalParamNEnd
    \cgalParamNBegin{max_octree_node_size}
      \cgalParamDescription{A node in the octree is only split if the contained number of primitives is larger and the maximal depth is not yet reached.}
      \cgalParamType{std::size_t}
      \cgalParamDefault{40}
    \cgalParamNEnd
  \cgalNamedParamsEnd


  \pre ! points.empty() and ! polygons.empty()

  */
  template<
    typename PointRange,
    typename PolygonRange,
    typename NamedParameters = parameters::Default_named_parameters>
  Kinetic_space_partition_3(
    const PointRange& points,
    const PolygonRange& polygons,
    const NamedParameters& np = CGAL::parameters::default_values()) :
    m_parameters(
      parameters::choose_parameter(parameters::get_parameter(np, internal_np::verbose), false),
      parameters::choose_parameter(parameters::get_parameter(np, internal_np::debug), false)), // use true here to export all steps
    m_input2regularized() {
    insert(points, polygons, np);
    initialize(np);
  }

  /*!
  \brief inserts non-coplanar polygons, requires `initialize()` afterwards to have effect.

  \tparam PointRange
  must be a model of `ConstRange` whose iterator type is `RandomAccessIterator` and whose value type is `GeomTraits::Point_3`.

  \tparam PolygonRange
  contains index ranges to form polygons by providing indices into PointRange

  \tparam NamedParameters
  a sequence of \ref bgl_namedparameters "Named Parameters"

  \param points
  an instance of `PointRange` with 3D points

  \param polygons
  a range of non-coplanar polygons defined by a range of indices into `points`

  \param np
  a sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

  \cgalNamedParamsBegin
     \cgalParamNBegin{point_map}
       \cgalParamDescription{a property map associating points to the elements of the `points`}
       \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type of the iterator of `PointRange` and whose value type is `GeomTraits::Point_3`}
       \cgalParamDefault{`CGAL::Identity_property_map<GeomTraits::Point_3>`}
     \cgalParamNEnd
  \cgalNamedParamsEnd
  */
  template<typename PointRange, typename PolygonRange, typename NamedParameters = parameters::Default_named_parameters>
  void insert(
    const PointRange& points,
    const PolygonRange& polygons,
    const NamedParameters& np = CGAL::parameters::default_values()) {

    using NP_helper = Point_set_processing_3_np_helper<PointRange, NamedParameters>;
    using PointMap = typename NP_helper::Point_map;

    PointMap point_map = NP_helper::get_point_map(np);

    To_exact to_exact;
    std::size_t offset = m_input2regularized.size();

    for (std::size_t p = 0; p < polygons.size();p++) {
      auto& poly = polygons[p];

      std::vector<Point_3> pts;
      pts.reserve(poly.size());
      for (auto it : poly)
        pts.push_back(get(point_map, *(points.begin() + it)));
      Plane_3 pl;
      Point_2 c;
      std::vector<Point_2> ch;
      process_input_polygon(pts, pl, c, ch);
      typename Intersection_kernel::Plane_3 exact_pl = to_exact(pl);

      // Check if there is already a coplanar polygon inserted
      bool skip = false;
      for (std::size_t i = 0; i < m_input_planes.size(); i++) {
        if (m_input_planes[i] == exact_pl) {
          std::cout << i << ". input polygon is coplanar to " << (p + offset) << ". input polygon" << std::endl;
          skip = true;
          break;
        }
      }

      if (skip)
        continue;

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

  /*!
  \brief initializes the kinetic partition of the bounding box.

  \tparam NamedParameters
  a sequence of \ref bgl_namedparameters "Named Parameters"

  \param np
  a sequence of \ref bgl_namedparameters "Named Parameters"
  among the ones listed below

  \cgalNamedParamsBegin
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
    \cgalParamNBegin{max_octree_depth}
      \cgalParamDescription{The maximal depth of the octree for subdividing the kinetic partition before initialization.}
      \cgalParamType{std::size_t}
      \cgalParamDefault{3}
    \cgalParamNEnd
    \cgalParamNBegin{max_octree_node_size}
      \cgalParamDescription{A node in the octree is only split if the contained number of primitives is larger and the maximal depth is not yet reached.}
      \cgalParamType{std::size_t}
      \cgalParamDefault{40}
    \cgalParamNEnd
  \cgalNamedParamsEnd

  \pre input data has been provided via `insert()`.
  */

  template<
    typename NamedParameters = parameters::Default_named_parameters>
  void initialize(
    const NamedParameters& np = CGAL::parameters::default_values()) {

    Timer timer;
    m_parameters.bbox_dilation_ratio = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::bbox_dilation_ratio), FT(11) / FT(10));
    m_parameters.reorient_bbox = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::reorient_bbox), false);
    m_parameters.max_octree_depth = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::max_octree_depth), 3);
    m_parameters.max_octree_node_size = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::max_octree_node_size), 40);

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
        KSP_3::internal::dump_polygon(m_input_polygons[i], std::to_string(i) + "-input_polygon");
    }

    split_octree();
    m_partitions.resize(m_partition_nodes.size());
    std::iota(m_partitions.begin(), m_partitions.end(), 0);

    for (std::size_t idx : m_partitions) {
      Sub_partition& partition = m_partition_nodes[idx];
      partition.index = idx;

      partition.m_data = std::make_shared<Data_structure>(m_parameters, std::to_string(idx) + "-");

      Initializer initializer(partition.clipped_polygons, partition.m_input_planes, *partition.m_data, m_parameters);
      initializer.initialize(partition.bbox, partition.input_polygons);
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

  \pre initialized partition and `k != 0`
  */
  void partition(std::size_t k) {
    FT a, b, c;
    partition(k, a, b, c);
  }

#ifndef DOXYGEN_RUNNING
  void partition(std::size_t k, FT& partition_time, FT& finalization_time, FT& conformal_time) {
    m_volumes.clear();
    Timer timer;
    timer.start();
    partition_time = 0;
    finalization_time = 0;
    conformal_time = 0;

    /*
        if (m_parameters.debug)
          if (boost::filesystem::is_directory("volumes/"))
            for (boost::filesystem::directory_iterator end_dir_it, it("volumes/"); it != end_dir_it; ++it)
              boost::filesystem::remove_all(it->path());*/

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
      Propagation propagation(*partition.m_data, m_parameters);
      std::size_t m_num_events = propagation.propagate(k);

      partition_time += timer.time();

      if (m_parameters.verbose) {
        std::cout << "* propagation finished" << std::endl;
        std::cout << "* number of events handled: " << m_num_events << std::endl;
      }

      if (m_parameters.verbose) {
        std::cout << std::endl << "--- FINALIZING PARTITION:" << std::endl;
      }

      // Finalization.

      for (std::size_t i = 0; i < partition.m_data->number_of_support_planes(); i++)
        if (!partition.m_data->support_plane(i).mesh().is_valid(true))
          std::cout << i << ". support has an invalid mesh!" << std::endl;

      for (std::size_t i = 6; i < partition.m_data->number_of_support_planes(); i++) {
        bool initial = false;
        typename Data_structure::Support_plane& sp = partition.m_data->support_plane(i);

        for (const auto& f : sp.mesh().faces())
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
    }

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

    std::map<typename Intersection_kernel::Point_3, Index> pts2idx;

    for (std::size_t i = 0; i < number_of_volumes(); i++) {
      std::vector<Index> f1;
      faces(i, std::back_inserter(f1));
      for (const Index& f : f1) {
        std::vector<Index>& face = m_partition_nodes[f.first].face2vertices[f.second];
        for (std::size_t j = 0; j < face.size(); j++) {
          auto it = pts2idx.emplace(m_partition_nodes[face[j].first].m_data->exact_vertices()[face[j].second], face[j]);
          if (!it.second)
            face[j] = it.first->second;
        }
      }
    }

    timer.stop();

 /*   if (m_parameters.debug) {

      if (boost::filesystem::is_directory("volumes/"))
        for (boost::filesystem::directory_iterator end_dir_it, it("volumes/"); it != end_dir_it; ++it)
          boost::filesystem::remove_all(it->path());

      KSP_3::dump_volumes_ksp(*this, "volumes/");
      for (std::size_t i = 1; i < m_volumes.size(); i++)
        if (m_volumes[i].first != m_volumes[i - 1].first)
          std::cout << i << " " << m_volumes[i - 1].first << std::endl;
      std::cout << m_volumes.size() << " " << m_volumes.back().first << std::endl;
    }*/

    timer.reset();
    timer.start();
    make_conformal(0);
    conformal_time = timer.time();

/*    if (m_parameters.debug) {

      if (boost::filesystem::is_directory("volumes_after/"))
        for (boost::filesystem::directory_iterator end_dir_it, it("volumes_after/"); it != end_dir_it; ++it)
          boost::filesystem::remove_all(it->path());
      KSP_3::dump_volumes_ksp(*this, "volumes_after/");
      for (std::size_t i = 1; i < m_volumes.size(); i++)
        if (m_volumes[i].first != m_volumes[i - 1].first)
          std::cout << i << " " << m_volumes[i - 1].first << std::endl;
      std::cout << m_volumes.size() << " " << m_volumes.back().first << std::endl;
    }*/

    if (m_parameters.verbose)
      check_tjunctions();

    // Clear unused data structures
    for (std::size_t i = 0; i < m_partitions.size(); i++) {
      m_partition_nodes[i].m_data->pface_neighbors().clear();
      m_partition_nodes[i].m_data->face_to_vertices().clear();
      m_partition_nodes[i].m_data->face_to_index().clear();
      m_partition_nodes[i].m_data->face_to_volumes().clear();
    }

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
  \brief returns the number of volumes created by the kinetic partition.

  \pre created partition
  */
  std::size_t number_of_volumes() const {
    return m_volumes.size();
  }

  /*!
  \brief provides the support planes of the partition derived from the input polygons

  @return
   vector of planes.

  \pre inserted polygons
  */
  const std::vector<typename Intersection_kernel::Plane_3> &input_planes() const {
    return m_input_planes;
  }

  /*!
   \brief exports the kinetic partition into a `Linear_cell_complex_for_combinatorial_map<3, 3>` using a model of `KineticLCCItems` as items, e.g., `Kinetic_space_partition_3::Linear_cell_complex_min_items`.

   Volume and face attributes defined in the model `KineticLCCItems` are filled. The volume index is in the range [0, number of volumes -1]

   \tparam LCC must be a model of `Linear_cell_complex_for_combinatorial_map<3, 3>` using a model of `KineticLCCItems`.

   \param lcc instance of LCC to be filled with the kinetic partition. Any data contained in `lcc` will be cleared before.

   \pre created partition
  */
  template<class LCC>
  void get_linear_cell_complex(LCC &lcc) const {
    lcc.clear();

    std::map<Index, std::size_t> mapped_vertices;
    std::map<typename Intersection_kernel::Point_3, std::size_t> mapped_points;
    std::vector<typename Intersection_kernel::Point_3> vtx;
    std::vector<Index> vtx_index;

    From_exact to_inexact;
    To_exact to_exact;

    std::vector<Index> faces_of_volume, vtx_of_face;
    std::vector<typename Intersection_kernel::Point_3> pts_of_face;

    for (std::size_t i = 0; i < number_of_volumes(); i++) {
      faces(i, std::back_inserter(faces_of_volume));

      for (const Index& f : faces_of_volume) {
        exact_vertices(f, std::back_inserter(pts_of_face), std::back_inserter(vtx_of_face));

        for (std::size_t j = 0; j < pts_of_face.size(); j++) {
          auto pit = mapped_points.emplace(pts_of_face[j], vtx.size());
          if (pit.second) {
            mapped_vertices[vtx_of_face[j]] = vtx.size();
            vtx.push_back(pts_of_face[j]);
            vtx_index.push_back(vtx_of_face[j]);
          }
          else mapped_vertices[vtx_of_face[j]] = pit.first->second;
        }

        pts_of_face.clear();
        vtx_of_face.clear();
      }
      faces_of_volume.clear();
    }

    CGAL::Linear_cell_complex_incremental_builder_3<LCC> ib(lcc);
    for (std::size_t i = 0; i < vtx.size(); i++)
      ib.add_vertex(vtx[i]);

    std::size_t num_faces = 0;
    std::size_t num_vols = 0;
    std::size_t num_vtx = 0;

    typename LCC::Dart_descriptor d;

    std::vector<bool> used_vertices(mapped_vertices.size(), false);
    std::vector<bool> added_volumes(number_of_volumes(), false);
    std::deque<std::size_t> queue;
    queue.push_back(0);
    while (!queue.empty()) {
      std::size_t v = queue.front();
      queue.pop_front();

      if (added_volumes[v])
        continue;

      if (!can_add_volume_to_lcc(v, added_volumes, mapped_vertices, used_vertices)) {
        queue.push_back(v);
        continue;
      }

      added_volumes[v] = true;

      ib.begin_surface();
      //std::cout << v << " inserting:";
      num_vols++;
      faces(v, std::back_inserter(faces_of_volume));

      typename Intersection_kernel::Point_3 centroid = to_exact(m_partition_nodes[m_volumes[v].first].m_data->volumes()[m_volumes[v].second].centroid);

      /*
            std::ofstream vout3(std::to_string(v) + ".xyz");
            vout3.precision(20);
            vout3 << " " << m_partition_nodes[m_volumes[v].first].m_data->volumes()[m_volumes[v].second].centroid << std::endl;
            vout3 << std::endl;
            vout3.close();*/

            // How to order faces accordingly?
            // First take faces of adjacent volumes and collect all added edges
            // Then pick from the remaining faces and take those which have already inserted edges
            // Repeat the last step until all are done.
      //       std::set<std::pair<std::size_t, std::size_t> > edges;
      //       for (std::size_t j=0;)
            // Try easy way and remove cells, I did not add after every loop?

      for (std::size_t j = 0; j < faces_of_volume.size(); j++) {
        vertex_indices(faces_of_volume[j], std::back_inserter(vtx_of_face));

        auto pair = neighbors(faces_of_volume[j]);

        if (pair.first != static_cast<int>(v) && !added_volumes[pair.first])
          queue.push_back(pair.first);
        if (pair.second != static_cast<int>(v) && pair.second >= 0 && !added_volumes[pair.second])
          queue.push_back(pair.second);

        //auto vertex_range = m_data.pvertices_of_pface(vol.pfaces[i]);
        ib.begin_facet();
        num_faces++;

        //std::cout << "(";

        //Sub_partition& p = m_partition_nodes[faces_of_volume[j].first];

        typename Intersection_kernel::Vector_3 norm;
        std::size_t i = 0;
        do {
          std::size_t n = (i + 1) % vtx_of_face.size();
          std::size_t nn = (n + 1) % vtx_of_face.size();
          norm = CGAL::cross_product(vtx[mapped_vertices[vtx_of_face[n]]] - vtx[mapped_vertices[vtx_of_face[i]]], vtx[mapped_vertices[vtx_of_face[nn]]] - vtx[mapped_vertices[vtx_of_face[n]]]);
          i++;
        } while (to_inexact(norm.squared_length()) == 0 && i < vtx_of_face.size());

        FT len = sqrt(to_inexact(norm.squared_length()));
        if (len != 0)
          len = 1.0 / len;
        norm = norm * to_exact(len);

        bool outwards_oriented = (vtx[mapped_vertices[vtx_of_face[0]]] - centroid) * norm < 0;
        //outward[std::make_pair(v, j)] = outwards_oriented;

        if (!outwards_oriented)
          std::reverse(vtx_of_face.begin(), vtx_of_face.end());

        /*
                auto p1 = edge_to_volface.emplace(std::make_pair(std::make_pair(mapped_vertices[vtx_of_face[0]], mapped_vertices[vtx_of_face[1]]), std::make_pair(v, j)));
                if (!p1.second) {
                  std::size_t first = mapped_vertices[vtx_of_face[0]];
                  std::size_t second = mapped_vertices[vtx_of_face[1]];
                  auto p = edge_to_volface[std::make_pair(first, second)];
                  auto o1 = outward[p];
                  auto o2 = outward[std::make_pair(v, j)];
                }

                for (std::size_t k = 1; k < vtx_of_face.size() - 1; k++) {
                  auto p = edge_to_volface.emplace(std::make_pair(std::make_pair(mapped_vertices[vtx_of_face[k]], mapped_vertices[vtx_of_face[k + 1]]), std::make_pair(v, j)));
                  if (!p.second) {
                    std::size_t first = mapped_vertices[vtx_of_face[k]];
                    std::size_t second = mapped_vertices[vtx_of_face[k + 1]];
                    auto p = edge_to_volface[std::make_pair(first, second)];
                    auto o1 = outward[p];
                    auto o2 = outward[std::make_pair(v, j)];
                  }
                }

                auto p2 = edge_to_volface.emplace(std::make_pair(std::make_pair(mapped_vertices[vtx_of_face.back()], mapped_vertices[vtx_of_face[0]]), std::make_pair(v, j)));
                if (!p2.second) {
                  std::size_t first = mapped_vertices[vtx_of_face.back()];
                  std::size_t second = mapped_vertices[vtx_of_face[0]];
                  auto p = edge_to_volface[std::make_pair(first, second)];
                  auto o1 = outward[p];
                  auto o2 = outward[std::make_pair(v, j)];
                }*/

        for (const Index& v : vtx_of_face) {
          ib.add_vertex_to_facet(static_cast<std::size_t>(mapped_vertices[v]));
          //std::cout << " " << mapped_vertices[v];
          if (!used_vertices[mapped_vertices[v]]) {
            used_vertices[mapped_vertices[v]] = true;
            num_vtx++;
          }
        }

        //std::cout << ")";
        auto face_dart = ib.end_facet(); // returns a dart to the face
        if (lcc.template attribute<2>(face_dart) == lcc.null_descriptor) {
          lcc.template set_attribute<2>(face_dart, lcc.template create_attribute<2>());
          // How to handle bbox planes that coincide with input polygons? Check support plane
          std::size_t sp = m_partition_nodes[faces_of_volume[j].first].m_data->face_to_support_plane()[faces_of_volume[j].second];

          // There are three different cases:
          // 1. face belongs to a plane from an input polygon
          // 2. face originates from octree splitting (and does not have an input plane)
          // 3. face lies on the bbox
          int ip = static_cast<int>(m_partition_nodes[faces_of_volume[j].first].m_data->support_plane(sp).data().actual_input_polygon);
          if (ip != -1)
            lcc.template info<2>(face_dart).input_polygon_index = static_cast<Face_support>(m_partition_nodes[faces_of_volume[j].first].input_polygons[ip]);
          else {
            // If there is no input polygon, I need to check whether it has two neighbors
            auto n = neighbors(faces_of_volume[j]);
            if (n.second >= 0)
              lcc.template info<2>(face_dart).input_polygon_index = static_cast<Face_support>(-7);
            else
              lcc.template info<2>(face_dart).input_polygon_index = static_cast<Face_support>(n.second);
          }
          lcc.template info<2>(face_dart).part_of_initial_polygon = m_partition_nodes[faces_of_volume[j].first].m_data->face_is_part_of_input_polygon()[faces_of_volume[j].second];
          lcc.template info<2>(face_dart).plane = m_partition_nodes[faces_of_volume[j].first].m_data->support_plane(m_partition_nodes[faces_of_volume[j].first].m_data->face_to_support_plane()[faces_of_volume[j].second]).exact_plane();
        }
        else {
          assert(lcc.template info<2>(face_dart).part_of_initial_polygon == m_partition_nodes[faces_of_volume[j].first].m_data->face_is_part_of_input_polygon()[faces_of_volume[j].second]);
        }

        vtx_of_face.clear();
      }

      d = ib.end_surface();

      auto ah = lcc.template create_attribute<3>();
      lcc.template set_attribute<3>(d, ah);
      lcc.template info<3>(d).barycenter = centroid;
      lcc.template info<3>(d).volume_id = v;

      faces_of_volume.clear();
    }

    // Todo: Remove check if all volumes were added
    for (std::size_t i = 0; i < added_volumes.size(); i++)
      if (!added_volumes[i])
        std::cout << "volume " << i << " has not been added" << std::endl;

    std::cout << "lcc #volumes: " << lcc.template one_dart_per_cell<3>().size() << " ksp #volumes: " << number_of_volumes() << std::endl;
    std::cout << "lcc #faces: " << lcc.template one_dart_per_cell<2>().size() << " ksp #faces: " << num_faces << std::endl;
    std::cout << "lcc #n-edges: " << lcc.template one_dart_per_cell<1>().size() << std::endl;
    std::cout << "lcc #vtx: " << lcc.template one_dart_per_cell<0>().size() << " ksp #vtx: " << vtx.size() << std::endl;

    // Verification
    // Check attributes in each dart
    for (auto& d : lcc.template one_dart_per_cell<0>()) {
      if (!lcc.is_dart_used(lcc.dart_descriptor(d))) {
        std::cout << "unused dart in 0" << std::endl;
      }
    }
    for (auto& d : lcc.template one_dart_per_cell<1>()) {
      if (!lcc.is_dart_used(lcc.dart_descriptor(d))) {
        std::cout << "unused dart in 1" << std::endl;
      }
    }
    for (auto& d : lcc.template one_dart_per_cell<2>()) {
      if (!lcc.is_dart_used(lcc.dart_descriptor(d))) {
        std::cout << "unused dart in 2" << std::endl;
      }
    }
    for (auto& d : lcc.template one_dart_per_cell<3>()) {
      if (!lcc.is_dart_used(lcc.dart_descriptor(d))) {
        std::cout << "unused dart in 3" << std::endl;
      }
    }

    lcc.display_characteristics(std::cout) << std::endl;

    if (!lcc.is_valid())
      std::cout << "LCC is not valid" << std::endl;
  }

  /// @}

private:
  struct Constraint_info {
    typename CDTplus::Constraint_id id_single, id_merged, id_overlay;
    std::size_t volume;
    Index vA, vB;
  };

  const Point_3& volume_centroid(std::size_t volume_index) const {
    assert(volume_index < m_volumes.size());
    auto p = m_volumes[volume_index];
    return m_partition_nodes[p.first].m_data->volumes()[p.second].centroid;
  }


  /*!
  \brief Face indices of the volume.

  \param volume_index
   index of the query volume.

  @return
   vector of face indices.

  \pre created partition
  */
  template<class OutputIterator>
  void faces(std::size_t volume_index, OutputIterator it) const {
    CGAL_assertion(m_volumes.size() > volume_index);
    auto p = m_volumes[volume_index];

    for (std::size_t i : m_partition_nodes[p.first].m_data->volumes()[p.second].faces)
      *it++ = std::make_pair(p.first, i);
  }


  /*!
  \brief Mapping of a vertex index to its position.

  @return
   vector of points.

    \pre created partition
  */
  const Point_3& vertex(const Index& vertex_index) const {
    return m_partition_nodes[vertex_index.first].m_data->vertices()[vertex_index.second];
  }

  /*!
  \brief Mapping of a vertex index to its exact position.

  @return
   vector of points.

    \pre created partition
  */
  const typename Intersection_kernel::Point_3& exact_vertex(const Index& vertex_index) const {
    return m_partition_nodes[vertex_index.first].m_data->exact_vertices()[vertex_index.second];
  }

  /*!
  \brief Vertices of a face.

  \param volume_index
   index of the query volume.

  @return
   vector of face indices.

  \pre created partition
  */
  template<class OutputIterator>
  void vertices(const Index& face_index, OutputIterator it) const {
    for (auto& p : m_partition_nodes[face_index.first].face2vertices[face_index.second])
      *it++ = m_partition_nodes[p.first].m_data->vertices()[p.second];
  }

  template<class OutputIterator>
  void vertex_indices(const Index& face_index, OutputIterator it) const {
    for (auto& p : m_partition_nodes[face_index.first].face2vertices[face_index.second])
      *it++ = p;
  }

  /*!
  \brief Vertices of a face.

  \param volume_index
   index of the query volume.

  @return
   vector of face indices.

  \pre created partition
  */
  template<class OutputIterator>
  void exact_vertices(const Index& face_index, OutputIterator it) const {

    for (auto& p : m_partition_nodes[face_index.first].face2vertices[face_index.second])
      *it++ = m_partition_nodes[p.first].m_data->exact_vertices()[p.second];
  }

  template<class OutputIterator, class IndexOutputIterator>
  void exact_vertices(const Index& face_index, OutputIterator pit, IndexOutputIterator iit) const {
    for (auto& p : m_partition_nodes[face_index.first].face2vertices[face_index.second]) {
      *iit++ = p;
      *pit++ = m_partition_nodes[p.first].m_data->exact_vertices()[p.second];
    }
  }


  template<class OutputIterator>
  void faces_of_input_polygon(const std::size_t polygon_index, OutputIterator it) const {
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
      const auto& f2sp = p.m_data->face_to_support_plane();

      for (std::size_t i = 0; i < f2sp.size(); i++) {
        if (f2sp[i] == sp_idx)
          *it++ = std::make_pair(idx, i);
      }
    }
  }

  const std::vector<std::vector<std::size_t> >& input_mapping() const {
    return m_regularized2input;
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

  \pre created partition
  */
  std::pair<int, int> neighbors(const Index &face_index) const {
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

  std::pair<int, int> make_canonical_pair(int i, int j) {
    if (i > j) return std::make_pair(j, i);
    return std::make_pair(i, j);
  }

  double build_cdt(CDTplus& cdt, std::vector<Index>& faces, std::vector<std::vector<Constraint_info> >& constraints, const typename Intersection_kernel::Plane_3& plane) {
    double area = 0;
    From_exact from_exact;

    cdt.clear();
    //keep track of constraints when inserting to iterate later
    constraints.resize(faces.size());

    //check orientation of faces so that they are ccw oriented
    std::vector<std::vector<Index> > pts_idx(faces.size());
    std::vector<std::vector<typename Intersection_kernel::Point_3> > pts(faces.size());
    for (std::size_t i = 0; i < faces.size(); ++i) {
      exact_vertices(faces[i], std::back_inserter(pts[i]), std::back_inserter(pts_idx[i]));
      constraints[i].resize(pts[i].size());

      CGAL::Orientation res = CGAL::COLLINEAR;
      bool pos = false;
      bool neg = false;

      for (std::size_t j = 0; j < pts[i].size(); j++) {
        std::size_t k = (j + 1) % pts[i].size();
        std::size_t l = (k + 1) % pts[i].size();

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
        vertices.push_back(cdt.insert(plane.to_2d(pts[f][v])));

        if (vertices.back()->info().idA2.first != static_cast<std::size_t>(-1) && vertices.back()->info().idA2 != pts_idx[f][v]) {
          std::cout << "build_cdt faces has non-unique vertices" << std::endl;
        }

        vertices.back()->info().idA2 = pts_idx[f][v];
        assert(pts_idx[f][v].first != static_cast<std::size_t>(-1));
        assert(pts_idx[f][v].second != static_cast<std::size_t>(-1));
        vertices.back()->info().adjacent.insert(faces[f]);
        vertices.back()->info().set_point(pts[f][v]);
        face2vtx[pts_idx[f][v]] = vertices.size() - 1;
        vtx2face[vertices.size() - 1] = pts_idx[f][v];
      }

    typedef std::set<std::pair<int, int> > Edges;
    Edges edges;

    // Iterating over each face and inserting each edge as a constraint.
    for (std::size_t i = 0; i < pts_idx.size(); ++i) {
      auto& v = pts_idx[i];
      for (std::size_t j = 0; j < v.size(); ++j) {
        int vj = static_cast<int>(face2vtx[v[j]]);
        int vjj = static_cast<int>(face2vtx[v[(j + 1) % v.size()]]);
        std::pair<Edges::iterator, bool> res = edges.insert(make_canonical_pair(vj, vjj));

        if (res.second) {
          constraints[i][j].id_single = cdt.insert_constraint(vertices[vj], vertices[vjj]);
          auto p = neighbors(faces[i]);
          if (p.second >= 0)
            std::cout << "p.second is positive" << std::endl;
          if (p.first < 0)
            std::cout << "p.first is negative" << std::endl;
          constraints[i][j].volume = p.first;
          constraints[i][j].vA = v[j];
          constraints[i][j].vB = v[(j + 1) % v.size()];
        }
      }
    }

    for (typename CDTplus::Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
      std::set<Index>& a(fit->vertex(0)->info().adjacent), & b(fit->vertex(1)->info().adjacent), & c(fit->vertex(2)->info().adjacent);

      std::set<Index> res, res2;
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

  void check_tjunctions() {
    std::map<Index, std::vector<Index> > vertex2neighbors;

    for (std::size_t v = 0; v < m_volumes.size(); v++) {
      auto &vp = m_volumes[v];
      for (const std::size_t f : m_partition_nodes[vp.first].m_data->volumes()[vp.second].faces) {
        auto& vtx = m_partition_nodes[vp.first].face2vertices[f];
        for (std::size_t i = 0; i < vtx.size(); i++) {
          vertex2neighbors[vtx[i]].push_back(vtx[(i + 1) % vtx.size()]);
          vertex2neighbors[vtx[i]].push_back(vtx[(i - 1 + vtx.size()) % vtx.size()]);
        }
      }
    }

    for (auto& p : vertex2neighbors) {
      typename Intersection_kernel::Point_3 a = m_partition_nodes[p.first.first].m_data->exact_vertices()[p.first.second];
      //Check pairwise collinear
      for (std::size_t i = 0; i < p.second.size(); i++) {
        typename Intersection_kernel::Point_3 b = m_partition_nodes[p.second[i].first].m_data->exact_vertices()[p.second[i].second];
        for (std::size_t j = i + 1; j < p.second.size(); j++) {
          if (p.second[i] == p.second[j])
            continue;
          typename Intersection_kernel::Point_3 c = m_partition_nodes[p.second[j].first].m_data->exact_vertices()[p.second[j].second];
          if (CGAL::collinear(a, b, c) && ((b - a) * (c - a) > 0)) {
            std::cout << "non-manifold v (" << p.first.first << ", " << p.first.second << ")" << std::endl;
            std::cout << " v (" << p.second[i].first << ", " << p.second[i].second << ")" << std::endl;
            std::cout << " v (" << p.second[j].first << ", " << p.second[j].second << ")" << std::endl;

            From_exact from_exact;

            std::ofstream vout2("a.xyz");
            vout2.precision(20);
            vout2 << from_exact(a) << std::endl;
            vout2.close();
            std::ofstream vout3("b.xyz");
            vout3.precision(20);
            vout3 << from_exact(b) << std::endl;
            vout3.close();
            std::ofstream vout4("c.xyz");
            vout4.precision(20);
            vout4 << from_exact(c) << std::endl;
            vout4.close();

            for (std::size_t v = 0; v < m_volumes.size(); v++) {
              auto& vp = m_volumes[v];
              for (const std::size_t f : m_partition_nodes[vp.first].m_data->volumes()[vp.second].faces) {
                auto& vtx = m_partition_nodes[vp.first].face2vertices[f];
                bool hasa = false, hasb = false, hasc = false;
                for (std::size_t k = 0; k < vtx.size(); k++) {
                  if (vtx[k] == p.first)
                    hasa = true;
                  if (vtx[k] == p.second[i])
                    hasb = true;
                  if (vtx[k] == p.second[j])
                    hasc = true;
                }

                if (hasa && (hasb || hasc)) {
                  const std::string vfilename = std::to_string(v) + " " + std::to_string(f) + "-non_manifold.polylines.txt";
                  std::ofstream vout(vfilename);
                  vout.precision(20);
                  vout << vtx.size() + 1;
                  for (const auto& v : vtx) {
                    vout << " " << from_exact(m_partition_nodes[v.first].m_data->exact_vertices()[v.second]);
                  }

                  vout << " " << from_exact(m_partition_nodes[vtx[0].first].m_data->exact_vertices()[vtx[0].second]);

                  vout << std::endl;
                  vout.close();
                }
              }
            }
            std::cout << std::endl;
          }
        }
      }
    }
  }

  void insert_map(const Index& a, const Index& b, std::map<Index, Index>& pm) const {
    if (a == b)
      return;

    Index target = b;
    auto it = pm.find(b);
    if (it != pm.end())
      target = it->second;
    pm[a] = target;
  }

  void build_cdt(CDTplus& cdt, std::vector<CDTplus>& partitions, std::vector<std::vector<std::vector<Constraint_info> > >& constraints,  const typename Intersection_kernel::Plane_3& plane) {
    if (partitions.size() == 0)
      return;

    for (std::size_t i = 0; i < partitions.size(); i++) {
      std::vector<Vertex_handle> vertices;
      vertices.reserve(6);

      for (std::size_t j = 0; j < constraints[i].size(); j++)
        for (std::size_t k = 0; k < constraints[i][j].size(); k++) {
          if (constraints[i][j][k].id_single == 0)
            continue;
          for (typename CDTplus::Vertices_in_constraint_iterator vi = partitions[i].vertices_in_constraint_begin(constraints[i][j][k].id_single); vi != partitions[i].vertices_in_constraint_end(constraints[i][j][k].id_single); vi++) {
            vertices.push_back(*vi);
          }

          // Insert constraints and replacing vertex handles in vector while copying data.
          VI tmp = vertices[0]->info();
          vertices[0] = cdt.insert(vertices[0]->point());
          vertices[0]->info() = tmp;

          tmp = vertices.back()->info();
          vertices.back() = cdt.insert(vertices.back()->point());
          vertices.back()->info() = tmp;

          constraints[i][j][k].id_merged = cdt.insert_constraint(vertices[0], vertices.back());

          vertices.clear();
        }
    }

    // Generate 3D points corresponding to the intersections
    std::size_t newpts = 0;
    for (typename CDTplus::Finite_vertices_iterator vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit) {
      if (!vit->info().input) {
        vit->info().point_3 = plane.to_3d(vit->point());
        vit->info().idA2 = vit->info().idB2 = Index(-1, -1);
        newpts++;
      }
    }

    for (typename CDTplus::Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
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
  }

  std::pair<double, double> overlay(CDTplus& cdtC, const CDTplus& cdtA, std::vector<std::vector<std::vector<Constraint_info> > >& constraints_a, const CDTplus& cdtB, std::vector<std::vector<std::vector<Constraint_info> > >& constraints_b, const typename Intersection_kernel::Plane_3& plane) {
    From_exact from_exact;
    //To_exact to_exact;
    std::pair<double, double> result;
    cdtC = cdtA;

    std::vector<Vertex_handle> vertices;
    vertices.reserve(2);

    for (std::size_t i = 0; i < constraints_a.size(); i++)
    for (std::size_t j = 0; j < constraints_a[i].size(); j++)
        for (std::size_t k = 0; k < constraints_a[i][j].size(); k++) {
          if (constraints_a[i][j][k].id_merged == 0) {
            if (constraints_a[i][j][k].id_single != 0)
              constraints_a[i][j][k].id_merged = constraints_a[i][j][k].id_single;
            else
              continue;
          }

          for (typename CDTplus::Vertices_in_constraint_iterator vi = cdtA.vertices_in_constraint_begin(constraints_a[i][j][k].id_merged); vi != cdtA.vertices_in_constraint_end(constraints_a[i][j][k].id_merged); vi++) {
            vertices.push_back(*vi);
          }

          // Insert constraints and replacing vertex handles in vector while copying data.
          VI tmp = vertices[0]->info();
          vertices[0] = cdtC.insert(vertices[0]->point());
          vertices[0]->info() = tmp;

          tmp = vertices.back()->info();
          vertices.back() = cdtC.insert(vertices.back()->point());
          vertices.back()->info() = tmp;

          constraints_a[i][j][k].id_overlay = cdtC.insert_constraint(vertices[0], vertices.back());

          vertices.clear();
        }

    for (std::size_t i = 0; i < constraints_b.size(); i++)
      for (std::size_t j = 0; j < constraints_b[i].size(); j++)
        for (std::size_t k = 0; k < constraints_b[i][j].size(); k++) {
          if (constraints_b[i][j][k].id_merged == 0) {
            if (constraints_b[i][j][k].id_single != 0)
              constraints_b[i][j][k].id_merged = constraints_b[i][j][k].id_single;
            else
              continue;
          }
          for (typename CDTplus::Vertices_in_constraint_iterator vi = cdtB.vertices_in_constraint_begin(constraints_b[i][j][k].id_merged); vi != cdtB.vertices_in_constraint_end(constraints_b[i][j][k].id_merged); vi++) {
            vertices.push_back(*vi);
          }

          // Insert constraints and replacing vertex handles in vector while copying data.
          VI tmp = vertices[0]->info();
          vertices[0] = cdtC.insert(vertices[0]->point());
          vertices[0]->info() = tmp;

          tmp = vertices.back()->info();
          vertices.back() = cdtC.insert(vertices.back()->point());
          vertices.back()->info() = tmp;

          constraints_b[i][j][k].id_overlay = cdtC.insert_constraint(vertices[0], vertices.back());

          vertices.clear();
        }

    std::size_t newpts = 0;
    // Generate 3D points corresponding to the intersections
    for (typename CDTplus::Finite_vertices_iterator vit = cdtC.finite_vertices_begin(); vit != cdtC.finite_vertices_end(); ++vit) {
      if (!vit->info().input) {
        vit->info().point_3 = plane.to_3d(vit->point());
        vit->info().idA2 = vit->info().idB2 = Index(-1, -1);
        newpts++;
      }
    }

    for (typename CDTplus::Finite_faces_iterator cit = cdtC.finite_faces_begin(); cit != cdtC.finite_faces_end(); ++cit) {
      double a = 0;
      cit->info().id2 = std::make_pair(-1, -1);

      typename Intersection_kernel::Point_2 p = CGAL::centroid(cit->vertex(0)->point(), cit->vertex(1)->point(), cit->vertex(2)->point());
      typename CDTplus::Face_handle fhA = cdtA.locate(p);

      if (cdtA.is_infinite(fhA))
        std::cout << "No face located in A: " << from_exact(plane.to_3d(p)) << std::endl;

      if (fhA->info().id2 != std::make_pair(std::size_t(-1), std::size_t(-1))) {
        cit->info().idA2 = fhA->info().id2;
        result.first += a;
      }
      else
        std::cout << "Face in A is missing ID " << from_exact(plane.to_3d(p)) << std::endl;

      typename CDTplus::Face_handle fhB = cdtB.locate(p);
      if (cdtB.is_infinite(fhB))
        std::cout << "No face located in B: " << from_exact(plane.to_3d(p)) << std::endl;

      if (fhB->info().id2 != std::make_pair(std::size_t(-1), std::size_t(-1))) {
        cit->info().idB2 = fhB->info().id2;
        result.second += a;
      }
      else
        std::cout << "Face in B is missing ID " << from_exact(plane.to_3d(p)) << std::endl;
    }

    return result;
  }

  void collect_faces(std::size_t partition_idx, std::size_t sp_idx, std::vector<Index>& faces, typename Intersection_kernel::Plane_3& plane) {
    Sub_partition& p = m_partition_nodes[partition_idx];

    plane = p.m_data->support_plane(sp_idx).data().exact_plane;

    const std::vector<std::size_t>& f2sp = p.m_data->face_to_support_plane();

    for (std::size_t i = 0; i < f2sp.size(); i++)
      if (f2sp[i] == sp_idx)
        faces.push_back(std::make_pair(partition_idx, i));
  }

  void collect_faces(Octree_node node, std::size_t dimension, bool lower, std::vector<Index>& faces, typename Intersection_kernel::Plane_3& plane) {
    // Collects boundary faces of node from its children.
    // dimension specifies the axis of the boundary face and lower determines if it is the lower of upper face of the cube on the axis.

    // Support plane indices:
    // xmin 4, xmax 2
    // ymin 1, ymax 3
    // zmin 0, zmax 5

    if (m_octree->is_leaf(node)) {
      // Mapping to partition is needed.
      std::size_t idx = m_node2partition[node];

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

  bool can_add_volume_to_lcc(std::size_t volume, const std::vector<bool>& added_volumes, const std::map<Index, std::size_t> &vtx2index, const std::vector<bool>& added_vertices) const {
    std::set<Index> vertices_of_volume;
    std::vector<Index> faces_of_volume;
    faces(volume, std::back_inserter(faces_of_volume));

    for (std::size_t i = 0; i < faces_of_volume.size(); i++) {
      std::vector<Index> vtx;
      auto n = neighbors(faces_of_volume[i]);
      int other = (n.first == static_cast<int>(volume)) ? n.second : n.first;
      if (other < 0 || !added_volumes[other])
        continue;
      vertex_indices(faces_of_volume[i], std::back_inserter(vtx));

      for (std::size_t j = 0; j < vtx.size(); j++)
        vertices_of_volume.insert(vtx[j]);
    }

    for (std::size_t i = 0; i < faces_of_volume.size(); i++) {
      auto n = neighbors(faces_of_volume[i]);
      int other = (n.first == static_cast<int>(volume)) ? n.second : n.first;
      if (other >= 0 && added_volumes[other])
        continue;
      std::vector<Index> vtx;
      vertex_indices(faces_of_volume[i], std::back_inserter(vtx));

      for (std::size_t j = 0; j < vtx.size(); j++) {
        auto it = vtx2index.find(vtx[j]);
        assert(it != vtx2index.end());
        if (vertices_of_volume.find(vtx[j]) == vertices_of_volume.end() && added_vertices[it->second])
          return false;
      }
    }

    return true;
  }

  CGAL::Aff_transformation_3<Kernel> get_obb2abb(const std::vector<std::vector<Point_3> > &polys) const {
    std::vector<Point_2> pts2d;
    std::size_t size = 0;

    for (std::size_t i = 0; i < polys.size(); i++)
      size += polys[i].size();
    pts2d.reserve(size);

    FT minz = (std::numeric_limits<FT>::max)(), maxz = -(std::numeric_limits<FT>::max)();
    for (std::size_t i = 0; i < polys.size(); i++)
      for (std::size_t j = 0; j < polys[i].size(); j++) {
      pts2d.push_back(Point_2(polys[i][j].x(), polys[i][j].y()));
      minz = (std::min<FT>)(minz, polys[i][j].z());
      maxz = (std::max<FT>)(maxz, polys[i][j].z());
    }

    std::vector<Point_2> ch;
    CGAL::convex_hull_2(pts2d.begin(), pts2d.end(), std::back_inserter(ch));

    std::vector<Point_2> bbox;
    bbox.reserve(4);
    CGAL::min_rectangle_2(ch.begin(), ch.end(), std::back_inserter(bbox));

    Vector_2 axis1 = bbox[0] - bbox[1];
    Vector_2 axis2 = bbox[1] - bbox[2];
    FT la = CGAL::sqrt(axis1.squared_length());
    axis1 = axis1 * (1.0 / la);
    FT lb = CGAL::sqrt(axis2.squared_length());
    axis2 = axis2 * (1.0 / lb);

    if (CGAL::abs(axis1.x()) < CGAL::abs(axis2.x())) {
      Vector_2 tmp = axis1;
      axis1 = axis2;
      axis2 = tmp;
    }

    if (0 > axis1.x())
      axis1 = -axis1;

    axis2 = Vector_2(-axis1.y(), axis1.x());

    FT rot[9];
    rot[0] = axis1.x();
    rot[1] = axis1.y();
    rot[2] = 0.0;
    rot[3] = -axis1.y();
    rot[4] = axis1.x();
    rot[5] = 0;
    rot[6] = 1.0;
    rot[7] = 0;
    rot[8] = 0;

    CGAL::Aff_transformation_3<Kernel> R(axis1.x(), axis1.y(), 0,
                                        -axis1.y(), axis1.x(), 0,
                                                 0,         0, 1.0);

    CGAL::Aff_transformation_3<Kernel> T(CGAL::TRANSLATION, - typename Kernel::Vector_3((bbox[0].x() + bbox[2].x()) * 0.5, (bbox[0].y() + bbox[2].y()) * 0.5, (maxz + minz) * 0.5));

    return R * T;
  }

  bool same_face(const Face_handle& a, const Face_handle& b) const {
    return (b->info().idA2 == a->info().idA2 && b->info().idB2 == a->info().idB2);
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
      if (vi.idA2.first < vi.idB2.first)
        vertices[i] = vi.idA2;
      else if (vi.idB2.first != static_cast<std::size_t>(-1))
        vertices[i] = vi.idB2;
      else {
        std::size_t vidx = m_partition_nodes[f.first].m_data->vertices().size();
        m_partition_nodes[f.first].m_data->vertices().push_back(from_exact(vi.point_3));
        m_partition_nodes[f.first].m_data->exact_vertices().push_back(vi.point_3);
        vertices[i] = vi.idA2 = std::make_pair(f.first, vidx);
      }
    }
  }

  void adapt_faces(const CDTplus& cdt) {
    std::set<Index> replacedA, replacedB;

    std::size_t extracted = 0;
    for (typename CDTplus::Face_handle fh : cdt.finite_face_handles()) {
      // when extracting each face, I only need to insert vertices, that don't exist on either side. Otherwise, I can just reference the vertex in the other partition.
      // using cit->info().id2.first as visited flag. -1 means not visited
      if (fh->info().id2.first != static_cast<std::size_t>(-1))
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
        if (cdt.is_infinite(fh->neighbor(static_cast<int>(i))) || !same_face(fh, fh->neighbor(static_cast<int>(i)))) {
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
          auto first = eit;

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

  std::pair<std::size_t, int> find_portal(std::size_t volume, std::size_t former, const Index& vA, const Index& vB, std::size_t& portal) const {
    portal = static_cast<std::size_t>(-7);
    auto vol = m_volumes[volume];
    std::vector<std::size_t>& faces = m_partition_nodes[vol.first].m_data->volumes()[vol.second].faces;

    for (std::size_t f = 0; f < faces.size(); f++) {
      auto n = neighbors(std::make_pair(vol.first, faces[f]));
      if (n.first == static_cast<int>(former) || n.second == static_cast<int>(former))
        continue;

      std::size_t idxA = static_cast<std::size_t>(-1);
      std::size_t numVtx = m_partition_nodes[vol.first].face2vertices[faces[f]].size();
      for (std::size_t v = 0; v < numVtx; v++)
        if (m_partition_nodes[vol.first].face2vertices[faces[f]][v] == vA) {
          idxA = v;
          break;
        }
      // If vertex wasn't found, skip to next face.
      if (idxA == static_cast<std::size_t>(-1))
        continue;

      std::size_t idxB = static_cast<std::size_t>(-1);
      int dir = 0;
      if (m_partition_nodes[vol.first].face2vertices[faces[f]][(idxA + 1) % numVtx] == vB) {
        dir = 1;
        idxB = (idxA + 1) % numVtx;
      }
      else if (m_partition_nodes[vol.first].face2vertices[faces[f]][(idxA + numVtx - 1) % numVtx] == vB) {
        dir = -1;
        idxB = (idxA + numVtx - 1) % numVtx;
      }

      // If only the first vertex was found, it is just an adjacent face.
      if (idxB == static_cast<std::size_t>(-1))
        continue;

      // Edge found
      // Save portal face for next volume.
      portal = f;

      return std::make_pair(idxA, dir);
    }
    return std::make_pair(-1, -1);
  }

  void adapt_internal_edges(const CDTplus& cdtC, const std::vector<Index> &faces_node, std::vector<std::vector<Constraint_info> >& c) {
    assert(faces_node.size() == c.size());

    std::size_t not_skipped = 0;

    for (std::size_t f = 0; f < c.size(); f++) {
      std::vector<Index> faces_of_volume;
      // The face index is probably no longer valid and the full face has been replaced by a smaller face using merged indices
      // Each constraint has a volume.
      // Constraints of the same volume are subsequent
      for (std::size_t e = 0; e < c[f].size(); e++) {
        auto id = c[f][e].id_single;
        if (id == 0)
          continue;

        id = (c[f][e].id_merged != 0) ? c[f][e].id_merged : id;
        id = (c[f][e].id_overlay != 0) ? c[f][e].id_overlay : id;

        int volume = static_cast<int>(c[f][e].volume);

        //auto it = (c[f][e].vA < c[f][e].vB) ? constraint2edge.find(std::make_pair(c[f][e].vA, c[f][e].vB)) : constraint2edge.find(std::make_pair(c[f][e].vB, c[f][e].vA));

        // Extract edge
        std::vector<Index> vertices_of_edge;
        for (typename CDTplus::Vertices_in_constraint_iterator vi = cdtC.vertices_in_constraint_begin(id); vi != cdtC.vertices_in_constraint_end(id); vi++) {
          if ((*vi)->info().idA2.first == static_cast<std::size_t>(-1))
            vertices_of_edge.push_back((*vi)->info().idB2);
          else vertices_of_edge.push_back((*vi)->info().idA2);
        }

        // Not necessary, as I am replacing vertices anyway?
        if (vertices_of_edge.size() == 2)
          continue;

        not_skipped++;

        // Check length of constraint
        // size 2 means it has not been split, thus there are no t-junctions.
        assert (vertices_of_edge.size() >= 2);

        faces_of_volume.clear();
        faces(volume, std::back_inserter(faces_of_volume));

        int starting_volume = volume;

        std::size_t idx, idx2;
        auto p = find_portal(volume, -7, c[f][e].vA, c[f][e].vB, idx);

        if (idx == static_cast<std::size_t>(-7)) {
          continue;
        }
        auto n = neighbors(faces_of_volume[idx]);
        int other = (n.first == volume) ? n.second : n.first;
        auto p2 = find_portal(volume, other, c[f][e].vA, c[f][e].vB, idx2);

        // For cdtA, there should be two portals and for cdtB only one
        // How to discard the traversing one?
        if (idx != static_cast<std::size_t>(-7)) {
          // Check if the portal idx is traversing.
          // The neighbors of a portal can be negative if it is not in the current face between the octree nodes.

          if (idx2 < static_cast<std::size_t>(-7) && m_volumes[volume].first != m_volumes[other].first) {
            idx = idx2;
            p = p2;
          }
        }
        else {
          idx = idx2;
          p = p2;
        }
        if (idx == static_cast<std::size_t>(-7))
          continue;

        std::vector<Index> tmp = m_partition_nodes[faces_of_volume[idx].first].face2vertices[faces_of_volume[idx].second];

        // Insert vertices in between
        if (p.second == 1)
          for (std::size_t i = 1; i < vertices_of_edge.size() - 1; i++)
            m_partition_nodes[faces_of_volume[idx].first].face2vertices[faces_of_volume[idx].second].insert(m_partition_nodes[faces_of_volume[idx].first].face2vertices[faces_of_volume[idx].second].begin() + p.first + i, vertices_of_edge[i]);
        else
          for (std::size_t i = 1; i < vertices_of_edge.size() - 1; i++)
            m_partition_nodes[faces_of_volume[idx].first].face2vertices[faces_of_volume[idx].second].insert(m_partition_nodes[faces_of_volume[idx].first].face2vertices[faces_of_volume[idx].second].begin() + p.first, vertices_of_edge[i]);


        n = neighbors(faces_of_volume[idx]);

        if (n.first != volume && n.second != volume)
          std::cout << "portal does not belong to volume" << std::endl;
        volume = (n.first == volume) ? n.second : n.first;
        int former = (idx == idx2) ? -1 : static_cast<int>(idx2);

        while (volume >= 0 && volume != starting_volume) { // What are the stopping conditions? There are probably several ones, e.g., arriving at the starting volume, not finding a portal face
          faces_of_volume.clear();
          faces(volume, std::back_inserter(faces_of_volume));

          auto p = find_portal(volume, former, c[f][e].vA, c[f][e].vB, idx);

          if (idx == static_cast<std::size_t>(-7))
            break;

          // Insert vertices in between
          if (p.second == 1)
            for (std::size_t i = 1; i < vertices_of_edge.size() - 1; i++)
              m_partition_nodes[faces_of_volume[idx].first].face2vertices[faces_of_volume[idx].second].insert(m_partition_nodes[faces_of_volume[idx].first].face2vertices[faces_of_volume[idx].second].begin() + p.first + i, vertices_of_edge[i]);
          else
            for (std::size_t i = 1; i < vertices_of_edge.size() - 1; i++)
              m_partition_nodes[faces_of_volume[idx].first].face2vertices[faces_of_volume[idx].second].insert(m_partition_nodes[faces_of_volume[idx].first].face2vertices[faces_of_volume[idx].second].begin() + p.first, vertices_of_edge[i]);


          // This is redundant to get next?
          auto n = neighbors(faces_of_volume[idx]);

          if (n.first != volume && n.second != volume)
            std::cout << "portal does not belong to volume" << std::endl;
          volume = (n.first == volume) ? n.second : n.first;

          former = volume;
        }
      }
    }
  }

  void make_conformal(std::vector<Index>& a, std::vector<Index>& b, typename Intersection_kernel::Plane_3& plane) {
    std::unordered_map<std::size_t, std::vector<Index> > a_sets, b_sets;
    for (const Index& i : a)
      a_sets[i.first].push_back(i);
    for (const Index& i : b)
      b_sets[i.first].push_back(i);

    // At first, one CDTplus is created for each partition node
    std::vector<CDTplus> a_cdts(a_sets.size()), b_cdts(b_sets.size());

    std::vector< std::vector<std::vector<Constraint_info> > > a_constraints;
    std::vector< std::vector<std::vector<Constraint_info> > > b_constraints;

    std::map<Index, Index> point_mapping;

    std::size_t idx = 0;
    a_constraints.resize(a_sets.size());

    std::set<std::size_t> partitions;
    for (auto& p : a_sets) {
      partitions.insert(p.first);
      build_cdt(a_cdts[idx], p.second, a_constraints[idx], plane);

      idx++;
    }

    idx = 0;
    b_constraints.resize(b_sets.size());
    for (auto& p : b_sets) {
      partitions.insert(p.first);
      build_cdt(b_cdts[idx], p.second, b_constraints[idx], plane);

      idx++;
    }

    CDTplus cdtA, cdtB, cdtC;
    build_cdt(cdtA, a_cdts, a_constraints, plane);

    build_cdt(cdtB, b_cdts, b_constraints, plane);

    overlay(cdtC, cdtA, a_constraints, cdtB, b_constraints, plane);

    adapt_faces(cdtC);

    idx = 0;
    for (auto& p : a_sets) {
      adapt_internal_edges(cdtC, p.second, a_constraints[idx]);
      idx++;
    }

    idx = 0;
    for (auto& p : b_sets) {
      adapt_internal_edges(cdtC, p.second, b_constraints[idx]);
      idx++;
    }
  }

  void make_conformal(Octree_node node)  {
    // pts2index maps exact points to their indices with one unique index.
    // index_map maps indices to unique indices. Used to map the points inside a partition to a unique index.

    // Nothing to do for a leaf node.
    if (m_octree->is_leaf(node))
      return;

    // Make childs conformal
    for (std::size_t i = 0; i < 8; i++)
      make_conformal(m_octree->child(node, i));

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

    }
  }

  void split_octree() {
    // Octree creation for sub partition
    std::size_t count = 0;
    for (const auto& p : m_input_polygons)
      count += p.size();

    m_points.clear();
    m_points.reserve(count);
    m_polygons.reserve(m_input_polygons.size());

    To_exact to_exact;
    From_exact from_exact;

    if (m_parameters.reorient_bbox) {

      m_transform = to_exact(get_obb2abb(m_input_polygons));

      for (const auto& p : m_input_polygons) {
        std::size_t idx = m_points.size();
        for (const Point_3& pt : p)
          m_points.push_back(from_exact(m_transform.transform(to_exact(pt))));

        m_polygons.push_back(std::vector<std::size_t>(p.size()));
        std::iota(m_polygons.back().begin(), m_polygons.back().end(), idx);
      }
    }
    else {
      m_transform = CGAL::Aff_transformation_3<Intersection_kernel>(CGAL::IDENTITY);

      for (const auto& p : m_input_polygons) {
        std::size_t idx = m_points.size();
        std::copy(p.begin(), p.end(), std::back_inserter(m_points));
        m_polygons.push_back(std::vector<std::size_t>(p.size()));
        std::iota(m_polygons.back().begin(), m_polygons.back().end(), idx);
      }
    }

    m_octree = std::make_unique<Octree>(CGAL::Orthtree_traits_polygons<Kernel>(m_points, m_polygons, m_parameters.bbox_dilation_ratio));
    m_octree->refine(m_parameters.max_octree_depth, m_parameters.max_octree_node_size);

    std::size_t leaf_count = 0;
    std::size_t max_count = 0;

    for (typename Octree::Node_index node : m_octree->traverse(CGAL::Orthtrees::Leaves_traversal<Octree>(*m_octree))) {
      if (m_octree->is_leaf(node))
        leaf_count++;
      else
        std::cout << "Leaves_traversal traverses non-leaves" << std::endl;
      max_count = (std::max<std::size_t>)(max_count, node);
    }

    m_partition_nodes.resize(leaf_count);

    m_node2partition.resize(max_count + 1, std::size_t(-1));

    std::size_t idx = 0;
    CGAL::Aff_transformation_3<Intersection_kernel> inv = m_transform.inverse();
    for (typename Octree::Node_index node : m_octree->traverse(CGAL::Orthtrees::Leaves_traversal<Octree>(*m_octree)))
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

        if (m_parameters.reorient_bbox)
          for (std::size_t i = 0; i < 8; i++)
            m_partition_nodes[idx].bbox[i] = inv.transform(m_partition_nodes[idx].bbox[i]);

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
            m_partition_nodes[idx].clipped_polygons[i][j] = from_exact(inv.transform(to_exact(polys[i].second[j])));
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

          KSP_3::internal::dump_polygons(m_partition_nodes[idx].clipped_polygons, std::to_string(idx) + "-polys.ply");
        }
        idx++;
      }

    std::cout << "input split into " << m_partition_nodes.size() << " partitions" << std::endl;
  }
};

} // namespace CGAL

#endif // CGAL_KINETIC_SPACE_PARTITION_3_H
