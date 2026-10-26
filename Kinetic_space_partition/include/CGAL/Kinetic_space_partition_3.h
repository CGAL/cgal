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

#include <algorithm>
#include <array>
#include <atomic>
#include <numeric>

// CGAL includes.
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Real_timer.h>
#include <CGAL/tags.h>
#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/partitioner.h>
#endif
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
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
#include <CGAL/KSP_3/Vertex_propagation.h>
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
    using Use_index = CGAL::Tag_true;
    using Index_type = std::uint32_t;

    struct Face_attribute {
      int input_polygon_index; // Non-negative numbers represent the index of the input polygon. Negative numbers correspond to the values defined in the enum `Kinetic_space_partition_3::Face_support`.
      int input_plane_index; // Non-negative numbers represent the index of the plane of the input polygon into `input_planes()`. Negative numbers like input_polygon_index.
      typename Intersection_kernel::Plane_3 plane;
      bool part_of_initial_polygon;
    };

    struct Volume_attribute {
      typename Intersection_kernel::Point_3 barycenter;
      std::size_t volume_id;
    };

    template<class LCC>
    struct Dart_wrapper {
      using Vertex_cell_attribute = CGAL::Cell_attribute_with_point< LCC, void>;
      using Face_cell_attribute = CGAL::Cell_attribute<LCC, Face_attribute>;
      using Volume_cell_attribute = CGAL::Cell_attribute<LCC, Volume_attribute>;

      using Attributes = std::tuple<Vertex_cell_attribute, void, Face_cell_attribute, Volume_cell_attribute>;
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
  using IkPoint_3 = typename Intersection_kernel::Point_3;

  using Index = std::pair<std::size_t, std::size_t>;

  using Data_structure = KSP_3::internal::Data_structure<Kernel, Intersection_kernel>;

  using IVertex = typename Data_structure::IVertex;
  using IEdge   = typename Data_structure::IEdge;

  using From_exact = typename CGAL::Cartesian_converter<Intersection_kernel, Kernel>;
  using To_exact = typename CGAL::Cartesian_converter<Kernel, Intersection_kernel>;

  using Initializer = KSP_3::internal::Initializer<Kernel, Intersection_kernel>;
  using Propagation = KSP_3::internal::Vertex_propagation<Kernel, Intersection_kernel>;
  using Finalizer   = KSP_3::internal::Finalizer<Kernel, Intersection_kernel>;

  using Timer        = CGAL::Real_timer;
  using Parameters = KSP::internal::Parameters_3<FT>;

  using Octree = CGAL::Orthtree<CGAL::Orthtree_traits_polygons<Intersection_kernel> >;
  using Octree_node = typename Octree::Node_index;
  //using Octree = CGAL::Orthtree<CGAL::Orthtree_traits_polygons<Kernel> >;
  //using Octree_node = typename Octree::Node_index;

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
    Sub_partition() {}
    std::shared_ptr<Data_structure> m_data;
    std::array<typename Intersection_kernel::Point_3, 8> bbox;
    std::vector<std::size_t> input_polygons;
    std::vector<std::vector<typename Intersection_kernel::Point_3> > clipped_polygons;
    std::vector<typename Intersection_kernel::Plane_3> m_input_planes;
    std::size_t index;

    // Results of the finalization. They are extracted from m_data after the Finalizer has run and m_data is released afterwards.
    // They are the only data accessed by make_conformal() and the access functions.
    struct Volume {
      std::vector<std::size_t> faces; // Indices into face2vertices.
      Point_3 centroid;
    };

    std::vector<std::pair<Index, Index> > face_neighbors;
    std::vector<std::vector<Index> > face2vertices;
    std::vector<Volume> volumes;
    std::vector<Point_3> vertices;
    std::vector<IkPoint_3> exact_vertices;
    std::vector<std::size_t> face2sp; // Support plane index for each face.
    std::vector<std::size_t> face_input; // Index into input_polygons for each face, std::size_t(-1) if none.
    std::vector<bool> face_part_of_input_polygon;
    std::vector<typename Intersection_kernel::Plane_3> sp_exact_planes; // Exact plane for each support plane.
    std::vector<std::size_t> sp_input_polygon; // Index into input_polygons for each support plane, std::size_t(-1) if none.
    std::map<std::size_t, std::size_t> input2sp; // Maps index into input_polygons onto support plane index.

    typename Octree::Node_index node;
  };

  Parameters m_parameters;
  CGAL::Aff_transformation_3<Intersection_kernel> m_transform;
  std::vector<Sub_partition> m_partition_nodes; // Tree of partitions.
  std::vector<std::size_t> m_partitions; // Contains the indices of the leaf nodes, the actual partitions to be calculated.
  std::vector<Point_3> m_points;
  std::vector<std::vector<std::size_t> > m_polygons;
  std::vector<std::vector<Point_3> > m_input_polygons;
  std::vector<typename Intersection_kernel::Plane_3> m_input_planes;

  // Polygons that are coplanar up to a tolerance share a plane.
  struct Coplanar_group {
    std::vector<Point_3> points; // Vertices of the first polygon of the group.
    Plane_3 first_plane; // Plane fitted to the first polygon of the group.
    Plane_3 plane; // Plane fitted to the vertices of the polygons inserted together with the first polygon of the group.
  };
  std::vector<Coplanar_group> m_coplanar_groups;
  std::vector<typename Intersection_kernel::Plane_3> m_planes; // Plane of each group in m_coplanar_groups.
  std::vector<std::size_t> m_input_plane_indices; // Index into m_planes for each inserted polygon, std::size_t(-1) for dropped degenerate polygons.
  std::vector<std::size_t> m_input2regularized; // Mapping from inserted polygons to the input polygons of the partition, std::size_t(-1) for dropped degenerate polygons.
  std::vector<std::vector<std::size_t> > m_regularized2input; // Mapping from partitioning planes to original input polygons.
  std::unique_ptr<Octree> m_octree;
  std::vector<std::size_t> m_node2partition;

  std::vector<Index> m_volumes;
  std::map<Index, std::size_t> m_index2volume;

  bool m_kinetic_data_initialized = false; // True if initialize() created the kinetic data structures of all octree leaves.


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
    const NamedParameters& np = CGAL::parameters::default_values())
    : m_parameters(
      parameters::choose_parameter(parameters::get_parameter(np, internal_np::verbose), false),
      parameters::choose_parameter(parameters::get_parameter(np, internal_np::debug), false)), // use true here to export all steps
    m_input2regularized() {}

  /*!
  \brief constructs a kinetic space partition object.

  \tparam PointRange
  must be a model of `ConstRange` whose iterator type is `RandomAccessIterator` and whose value type is Point_3.

  \tparam PolygonRange
  contains index ranges to form polygons by providing indices into PointRange

  \tparam NamedParameters
  a sequence of \ref bgl_namedparameters "Named Parameters"

  \param points
  an instance of `PointRange` with 3D points and corresponding 3D normal vectors

  \param polygons
  a range of polygons defined by a range of indices into `points`

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
      \cgalParamDescription{Write timing and internal information to `std::cout`.}
      \cgalParamType{bool}
      \cgalParamDefault{false}
    \cgalParamNEnd
    \cgalParamNBegin{reorient_bbox}
      \cgalParamDescription{Use the oriented bounding box instead of the axis-aligned bounding box. While the z direction is maintained, the x axis is aligned with the largest variation in the horizontal plane.}
      \cgalParamType{bool}
      \cgalParamDefault{false}
    \cgalParamNEnd
    \cgalParamNBegin{bbox_dilation_ratio}
      \cgalParamDescription{Factor for extension of the bounding box of the input data to be used for the partition. Must be larger than 1.0.}
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
    const NamedParameters& np = CGAL::parameters::default_values())
    : m_parameters(
      parameters::choose_parameter(parameters::get_parameter(np, internal_np::verbose), false),
      parameters::choose_parameter(parameters::get_parameter(np, internal_np::debug), false)), // use true here to export all steps
    m_input2regularized() {
    insert(points, polygons, np);
  }

  /*!
  \brief inserts polygons, requires `initialize()` or `init_and_partition()` afterwards to have effect.

  A plane is fitted to the vertices of each polygon and the polygon is replaced by the convex hull of its vertices projected onto the plane.
  Polygons that are coplanar up to a small tolerance relative to their extent get the same plane, which is fitted to the vertices of all of them.
  This avoids very thin volumes between nearly coplanar polygons.
  Degenerate polygons, i.e., polygons with less than three vertices or with collinear vertices, are dropped.

  Coplanar polygons share a support plane in the partition and expand separately. Pieces of several coplanar polygons that lie in the same face of
  the arrangement of the intersection lines within the support plane at the beginning are merged into their convex hull. A face of the partition
  originating from several coplanar polygons refers to one of them as its input polygon, see `get_linear_cell_complex()`.

  \tparam PointRange
  must be a model of `ConstRange` whose iterator type is `RandomAccessIterator` and whose value type is `GeomTraits::Point_3`.

  \tparam PolygonRange
  contains index ranges to form polygons by providing indices into PointRange

  \tparam NamedParameters
  a sequence of \ref bgl_namedparameters "Named Parameters"

  \param points
  an instance of `PointRange` with 3D points

  \param polygons
  a range of polygons defined by a range of indices into `points`

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

    static_assert(std::is_same_v<typename NP_helper::Geom_traits, Kernel>);

    PointMap point_map = NP_helper::get_point_map(np);

    std::size_t offset = m_input2regularized.size();

    std::vector<std::vector<Point_3> > polygon_points(polygons.size());
    for (std::size_t p = 0; p < polygons.size(); p++) {
      polygon_points[p].reserve(polygons[p].size());
      for (auto it : polygons[p])
        polygon_points[p].push_back(get(point_map, *(points.begin() + it)));
    }

    // Polygons that are coplanar up to a tolerance are grouped. The polygons of a group share the plane fitted to the vertices
    // of all of them. Nearly coplanar polygons would otherwise create very thin volumes.
    // A polygon that is coplanar to a group of a former call of insert() uses the plane of that group.
    const std::size_t first_new_group = m_coplanar_groups.size();
    std::vector<std::size_t> group(polygons.size(), std::size_t(-1));
    for (std::size_t p = 0; p < polygons.size(); p++) {
      std::array<IkPoint_3, 3> span;
      // Degenerate polygons, i.e., with less than three vertices or collinear vertices, are dropped.
      if (!spanning_points(polygon_points[p], span))
        continue;

      Plane_3 pl;
      CGAL::linear_least_squares_fitting_3(polygon_points[p].begin(), polygon_points[p].end(), pl, CGAL::Dimension_tag<0>());

      for (std::size_t g = 0; g < m_coplanar_groups.size() && group[p] == std::size_t(-1); g++)
        if (nearly_coplanar(polygon_points[p], pl, m_coplanar_groups[g].points, m_coplanar_groups[g].first_plane))
          group[p] = g;

      if (group[p] == std::size_t(-1)) {
        group[p] = m_coplanar_groups.size();
        m_coplanar_groups.push_back(Coplanar_group{ polygon_points[p], pl, pl });
      }
    }

    std::vector<std::vector<Point_3> > group_points(m_coplanar_groups.size() - first_new_group);
    for (std::size_t p = 0; p < polygons.size(); p++)
      if (group[p] != std::size_t(-1) && group[p] >= first_new_group)
        std::copy(polygon_points[p].begin(), polygon_points[p].end(), std::back_inserter(group_points[group[p] - first_new_group]));

    To_exact to_exact;
    for (std::size_t g = 0; g < group_points.size(); g++) {
      CGAL::linear_least_squares_fitting_3(group_points[g].begin(), group_points[g].end(), m_coplanar_groups[first_new_group + g].plane, CGAL::Dimension_tag<0>());
      m_planes.push_back(to_exact(m_coplanar_groups[first_new_group + g].plane));
    }

    for (std::size_t p = 0; p < polygons.size(); p++) {
      if (group[p] == std::size_t(-1)) {
        m_input2regularized.push_back(std::size_t(-1));
        m_input_plane_indices.push_back(std::size_t(-1));
        continue;
      }

      const Plane_3& pl = m_coplanar_groups[group[p]].plane;

      std::vector<Point_2> ch;
      process_input_polygon(polygon_points[p], pl, ch);

      // Coplanar input polygons are inserted separately. They share a support plane in the kinetic partition.
      m_input2regularized.push_back(m_input_planes.size());
      m_regularized2input.push_back(std::vector<std::size_t>());
      m_regularized2input.back().push_back(p + offset);
      m_input_planes.push_back(m_planes[group[p]]);
      m_input_plane_indices.push_back(group[p]);
      m_input_polygons.push_back(std::vector<Point_3>(ch.size()));

      for (std::size_t i = 0; i < ch.size(); i++)
        m_input_polygons.back()[i] = pl.to_3d(ch[i]);
    }
  }

  /*!
  \brief initializes the kinetic partition.

  The input polygons are separated using an octree. Afterwards, the kinetic data structure of each octree leaf is initialized.
  The leaves are initialized independently and, if `ConcurrencyTag` is `Parallel_tag`, in parallel.

  The initialized kinetic data structures are kept. Hence, `partition()` can be called several times afterwards, e.g., with different values of `k`, without initializing again.

  \note The initialized kinetic data structures of all octree leaves exist at the same time. If only a single partition is needed,
  `init_and_partition()` has a considerably lower memory consumption.

  \tparam ConcurrencyTag enables sequential versus parallel algorithm. Possible values are `Sequential_tag`, `Parallel_tag`, and `Parallel_if_available_tag`.

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
      \cgalParamDescription{Factor for extension of the bounding box of the input data to be used for the partition. Must be larger than 1.0.}
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
  template<typename ConcurrencyTag = Parallel_if_available_tag,
    typename NamedParameters = parameters::Default_named_parameters>
  void initialize(
    const NamedParameters& np = CGAL::parameters::default_values()) {
#ifndef CGAL_LINKED_WITH_TBB
    static_assert (!std::is_convertible<ConcurrencyTag, Parallel_tag>::value,
                   "Parallel_tag is enabled but TBB is unavailable.");
#endif

    release_kinetic_data();
    clear_partition();

    if (!create_subpartitions(np))
      return;

    FT initialization_time = 0;
    if (!initialize_subpartitions<ConcurrencyTag>(initialization_time)) {
      release_kinetic_data();
      return;
    }

    m_kinetic_data_initialized = true;

    if (m_parameters.verbose)
      std::cout << "* initialization time of kinetic data structures: " << initialization_time << std::endl;
  }

  /*!
  \brief propagates the kinetic polygons in the initialized partition.

  The kinetic partitions of the octree leaves are computed independently and, if `ConcurrencyTag` is `Parallel_tag`, in parallel.
  The kinetic data structures initialized by `initialize()` are restored after the propagation of each leaf. Hence, `partition()` can be called several times with different values of `k`.
  The result of a previous call is replaced.

  \note In the parallel version, the state of several kinetic data structures is copied at the same time, one per running thread.
  Thus, the memory consumption is higher than in the sequential version.

  \tparam ConcurrencyTag enables sequential versus parallel algorithm. Possible values are `Sequential_tag`, `Parallel_tag`, and `Parallel_if_available_tag`.

  \param k
   maximum number of allowed intersections for each vertex of the polygon before its expansion stops.

  \pre `initialize()` has been called and `k != 0`
  */
  template<typename ConcurrencyTag = Parallel_if_available_tag>
  void partition(std::size_t k) {
    FT a, b, c;
    partition<ConcurrencyTag>(k, a, b, c);
  }

#ifndef DOXYGEN_RUNNING
  template<typename ConcurrencyTag = Parallel_if_available_tag>
  void partition(std::size_t k, FT& partition_time, FT& finalization_time, FT& conformal_time) {
#ifndef CGAL_LINKED_WITH_TBB
    static_assert (!std::is_convertible<ConcurrencyTag, Parallel_tag>::value,
                   "Parallel_tag is enabled but TBB is unavailable.");
#endif

    partition_time = 0;
    finalization_time = 0;
    conformal_time = 0;

    if (!m_kinetic_data_initialized) {
      std::cout << "Warning: initialize() has to be called before partition()!" << std::endl;
      return;
    }

    clear_partition();

    FT initialization_time = 0;
    if (!compute_subpartitions<ConcurrencyTag>(k, false, initialization_time, partition_time, finalization_time))
      return;

    merge_subpartitions(conformal_time);
  }
#endif

  /*!
  \brief initializes the kinetic partition and propagates the kinetic polygons.

  The input polygons are separated using an octree. Afterwards, the kinetic partition of each octree leaf is computed.
  Each leaf are computed independently and, if `ConcurrencyTag` is `Parallel_tag`, in parallel.

  In contrast to calling `initialize()` and `partition()`, the kinetic data structure of each leaf is created right before its propagation and released after its finalization.
  Hence, the memory consumption is considerably lower. To compute another partition, e.g., with a different `k`, `init_and_partition()` has to be called again.

  \note In the parallel version, several kinetic data structures exist at the same time, one per running thread.
  Thus, the memory consumption is higher than in the sequential version.

  \tparam ConcurrencyTag enables sequential versus parallel algorithm. Possible values are `Sequential_tag`, `Parallel_tag`, and `Parallel_if_available_tag`.

  \tparam NamedParameters
  a sequence of \ref bgl_namedparameters "Named Parameters"

  \param k
   maximum number of allowed intersections for each vertex of the polygon before its expansion stops.

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
      \cgalParamDescription{Factor for extension of the bounding box of the input data to be used for the partition. Must be larger than 1.0.}
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

  \pre input data has been provided via `insert()` and `k != 0`.
  */
  template<typename ConcurrencyTag = Parallel_if_available_tag,
    typename NamedParameters = parameters::Default_named_parameters>
  void init_and_partition(std::size_t k,
    const NamedParameters& np = CGAL::parameters::default_values()) {
    FT a, b, c;
    init_and_partition<ConcurrencyTag>(k, a, b, c, np);
  }

#ifndef DOXYGEN_RUNNING
  template<typename ConcurrencyTag = Parallel_if_available_tag,
    typename NamedParameters = parameters::Default_named_parameters>
  void init_and_partition(std::size_t k, FT& partition_time, FT& finalization_time, FT& conformal_time,
    const NamedParameters& np = CGAL::parameters::default_values()) {
#ifndef CGAL_LINKED_WITH_TBB
    static_assert (!std::is_convertible<ConcurrencyTag, Parallel_tag>::value,
                   "Parallel_tag is enabled but TBB is unavailable.");
#endif

    partition_time = 0;
    finalization_time = 0;
    conformal_time = 0;

    release_kinetic_data();
    clear_partition();

    if (!create_subpartitions(np))
      return;

    FT initialization_time = 0;
    CGAL::Real_timer t1;
    t1.start();
    if (!compute_subpartitions<ConcurrencyTag>(k, true, initialization_time, partition_time, finalization_time))
      return;
    t1.stop();
    std::cout << "p: " << t1.time() << "\n";
    if (m_parameters.verbose) {
      std::cout << "* initialization time of kinetic data structures: " << initialization_time << "\n";
      std::cout << "* partition time: " << partition_time << "\n";
      std::cout << "* finalization time: " << finalization_time << "\n";
    }

    t1.reset();
    t1.start();
    merge_subpartitions(conformal_time);
    t1.stop();
    std::cout << "c: " << t1.time() << std::endl;

    if (m_parameters.verbose)
      std::cout << "* conformal time: " << conformal_time << std::endl;
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
  \brief provides the planes of the input polygons.

  Polygons that are coplanar up to a small tolerance share a plane. Thus, each plane occurs only once.
  The planes are indexed by `input_plane_indices()` and by the input plane index stored in the face attributes of the linear cell complex, see `get_linear_cell_complex()`.

  @return
   vector of planes.

  \pre inserted polygons
  */
  const std::vector<typename Intersection_kernel::Plane_3> &input_planes() const {
    return m_planes;
  }

  /*!
  \brief provides the index of the plane of each input polygon into `input_planes()`.

  Polygons that are coplanar up to a small tolerance share a plane and thus have the same index.
  Degenerate polygons, i.e., polygons with less than three vertices or with collinear vertices, are dropped when inserted. Their index is `std::size_t(-1)`.

  @return
   vector of plane indices, one for each inserted polygon.

  \pre inserted polygons
  */
  const std::vector<std::size_t>& input_plane_indices() const {
    return m_input_plane_indices;
  }

  /*!
   \brief exports the kinetic partition into a `Linear_cell_complex_for_combinatorial_map<3, 3>` using a model of `KineticLCCItems` as items, e.g., `Kinetic_linear_cell_complex_min_items`.

   Volume and face attributes defined in the model `KineticLCCItems` are filled. The volume index is in the range [0, number of volumes -1]
   The input polygon index of a face refers to the input polygon the face originates from. If the face originates from several coplanar input polygons, one of them is provided.

   \tparam LCC must be a model of `Linear_cell_complex_for_combinatorial_map<3, 3>` using a model of `KineticLCCItems`.

   \param lcc instance of LCC to be filled with the kinetic partition. Any data contained in `lcc` will be cleared before.

   \pre created partition
  */
  template<class LCC>
  void get_linear_cell_complex(LCC &lcc) const {
    lcc.clear();

    std::map<Index, std::size_t> mapped_vertices;
    std::map<IkPoint_3, std::size_t> mapped_points;
    std::vector<IkPoint_3> vtx;
    std::vector<Index> vtx_index;

    using LCC_kernel = typename CGAL::Kernel_traits<typename LCC::Point>::Kernel;

    using To_lcc = CGAL::Cartesian_converter<Intersection_kernel, LCC_kernel>;

    To_lcc to_lcc;

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
      ib.add_vertex(to_lcc(vtx[i]));

    std::size_t num_faces = 0;

    typename LCC::Dart_descriptor d;

    std::vector<bool> used_vertices(mapped_vertices.size(), false);
    std::vector<bool> added_volumes(number_of_volumes(), false);

    std::deque<std::size_t> queue;
    if (number_of_volumes() > 0)
      queue.push_back(0);

    bool first = true;
    while (!queue.empty()) {
      std::size_t v = queue.front();
      queue.pop_front();

      if (added_volumes[v])
        continue;

      added_volumes[v] = true;

      ib.begin_surface();

      faces(v, std::back_inserter(faces_of_volume));

      IkPoint_3 centroid = to_exact(m_partition_nodes[m_volumes[v].first].volumes[m_volumes[v].second].centroid);

      std::vector<bool> face_added(faces_of_volume.size(), false);
      std::vector<std::vector<Index> > vtx_of_faces(faces_of_volume.size());

      // All vertices of the volume, used to orient its faces exactly.
      std::vector<std::size_t> vtx_of_volume;
      for (std::size_t i = 0; i < faces_of_volume.size(); i++) {
        vertex_indices(faces_of_volume[i], std::back_inserter(vtx_of_faces[i]));
        for (const Index& vi : vtx_of_faces[i])
          vtx_of_volume.push_back(mapped_vertices[vi]);
      }

      std::queue<std::size_t> face_queue;
      for (std::size_t i = 0;i<faces_of_volume.size();i++)
        face_queue.push(i);

      while (!face_queue.empty()) {
        std::size_t j = face_queue.front();
        face_queue.pop();

        if (vtx_of_faces[j].empty())
          vertex_indices(faces_of_volume[j], std::back_inserter(vtx_of_faces[j]));

        if (first) {
          add_face_to_lcc(lcc, ib, faces_of_volume[j], vtx_of_faces[j], vtx_of_volume, mapped_vertices, vtx, used_vertices);

          // Add neighboring faces to queue
          auto pair = neighbors(faces_of_volume[j]);

          if (pair.first != static_cast<int>(v) && !added_volumes[pair.first])
            queue.push_back(pair.first);
          if (pair.second != static_cast<int>(v) && pair.second >= 0 && !added_volumes[pair.second])
            queue.push_back(pair.second);

          face_added[j] = true;
          num_faces++;
          first = false;
        }
        else {
          if (face_queue.size() <= 3 || can_add_face_to_lcc(vtx_of_faces[j], mapped_vertices, used_vertices)) {
            add_face_to_lcc(lcc, ib, faces_of_volume[j], vtx_of_faces[j], vtx_of_volume, mapped_vertices, vtx, used_vertices);

            // Add neighboring faces to queue
            auto pair = neighbors(faces_of_volume[j]);

            if (pair.first != static_cast<int>(v) && !added_volumes[pair.first])
              queue.push_back(pair.first);
            if (pair.second != static_cast<int>(v) && pair.second >= 0 && !added_volumes[pair.second])
              queue.push_back(pair.second);

            face_added[j] = true;
            num_faces++;
            first = false;
          }
          else
            face_queue.push(j);
        }
      }

      d = ib.end_surface();

      auto ah = lcc.template create_attribute<3>();
      lcc.template set_attribute<3>(d, ah);
      lcc.template info<3>(d).barycenter = centroid;
      lcc.template info<3>(d).volume_id = v;

      faces_of_volume.clear();
    }

    // Check whether all volumes were added
    CGAL_assertion_code(
    for (std::size_t i = 0; i < added_volumes.size(); i++)
      if (!added_volumes[i])
        std::cout << "volume " << i << " has not been added" << std::endl;);

    if (m_parameters.verbose) {
      std::cout << "lcc #volumes: " << lcc.template one_dart_per_cell<3>().size() << " ksp #volumes: " << number_of_volumes() << std::endl;
      std::cout << "lcc #faces: " << lcc.template one_dart_per_cell<2>().size() << " ksp #faces: " << num_faces << std::endl;
      std::cout << "lcc #n-edges: " << lcc.template one_dart_per_cell<1>().size() << std::endl;
      std::cout << "lcc #vtx: " << lcc.template one_dart_per_cell<0>().size() << " ksp #vtx: " << vtx.size() << std::endl;
    }

    // Verification
    CGAL_assertion_code(
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
    });

    if (m_parameters.verbose)
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


  /*!
  \brief Face indices of the volume.

  \tparam OutputIterator output iterator accepting Index

  \param volume_index
   index of the query volume.

  \param it
  output iterator where faces of the volume are put

  \pre created partition
  */
  template<class OutputIterator>
  void faces(std::size_t volume_index, OutputIterator it) const {
    CGAL_assertion(m_volumes.size() > volume_index);
    auto p = m_volumes[volume_index];

    for (std::size_t i : m_partition_nodes[p.first].volumes[p.second].faces)
      *it++ = std::make_pair(p.first, i);
  }


  /*!
  \brief Mapping of a vertex index to its position.

  @return point of the vertex

    \pre created partition
  */
  const Point_3& vertex(const Index& vertex_index) const {
    return m_partition_nodes[vertex_index.first].vertices[vertex_index.second];
  }

  /*!
  \brief Vertices of a face.

  \tparam OutputIterator output iterator accepting `Point_3`

  \param face_index
   index of the query face.

  \param it output iterator where points of the face are put

  \pre created partition
  */
  template<class OutputIterator>
  void vertices(const Index& face_index, OutputIterator it) const {
    for (auto& p : m_partition_nodes[face_index.first].face2vertices[face_index.second])
      *it++ = m_partition_nodes[p.first].vertices[p.second];
  }

  template<class OutputIterator>
  void vertex_indices(const Index& face_index, OutputIterator it) const {
    for (auto& p : m_partition_nodes[face_index.first].face2vertices[face_index.second])
      *it++ = p;
  }

  /*!
  \brief Vertices of a face.

  \tparam OutputIterator output iterator accepting `EPECK::Point_3`

  \param face_index
   index of the query face.

  \param it output iterator where points of the face are put

  \pre created partition
  */
  template<class OutputIterator>
  void exact_vertices(const Index& face_index, OutputIterator it) const {

    for (auto& p : m_partition_nodes[face_index.first].face2vertices[face_index.second])
      *it++ = m_partition_nodes[p.first].exact_vertices[p.second];
  }

  template<class OutputIterator, class IndexOutputIterator>
  void exact_vertices(const Index& face_index, OutputIterator pit, IndexOutputIterator iit) const {
    for (auto& p : m_partition_nodes[face_index.first].face2vertices[face_index.second]) {
      *iit++ = p;
      *pit++ = m_partition_nodes[p.first].exact_vertices[p.second];
    }
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
      CGAL_assertion(it != m_index2volume.end());
      return std::pair<int, int>(static_cast<int>(it->second), static_cast<int>(p.second.second));
    }
    else {
      auto it1 = m_index2volume.find(p.first);
      CGAL_assertion(it1 != m_index2volume.end());
      auto it2 = m_index2volume.find(p.second);
      CGAL_assertion(it2 != m_index2volume.end());
      return std::pair<int, int>(static_cast<int>(it1->second), static_cast<int>(it2->second));
    }
  }

  // Computes the convex hull of the polygon projected onto the plane pl.
  void process_input_polygon(const std::vector<Point_3>& poly, const Plane_3& pl, std::vector<Point_2>& ch) const {
    std::vector<Point_2> pts2d(poly.size());
    for (std::size_t i = 0; i < poly.size(); i++)
      pts2d[i] = pl.to_2d(poly[i]);

    ch.clear();
    CGAL::convex_hull_2(pts2d.begin(), pts2d.end(), std::back_inserter(ch));
  }

  // Provides three points of the polygon that are not collinear. The points are chosen to span a large triangle.
  bool spanning_points(const std::vector<Point_3>& poly, std::array<IkPoint_3, 3>& span) const {
    if (poly.size() < 3)
      return false;

    std::size_t a = 0, b = 0, c = 0;
    FT max_dist = -1, max_area = -1;
    for (std::size_t i = 1; i < poly.size(); i++) {
      const FT d = CGAL::squared_distance(poly[a], poly[i]);
      if (d > max_dist) {
        max_dist = d;
        b = i;
      }
    }
    for (std::size_t i = 0; i < poly.size(); i++) {
      const FT area = CGAL::cross_product(poly[b] - poly[a], poly[i] - poly[a]).squared_length();
      if (area > max_area) {
        max_area = area;
        c = i;
      }
    }

    To_exact to_exact;
    span = { to_exact(poly[a]), to_exact(poly[b]), to_exact(poly[c]) };
    return !CGAL::collinear(span[0], span[1], span[2]);
  }

  // Checks whether two polygons are coplanar up to a tolerance relative to their extent, i.e.,
  // whether each vertex of a polygon is close to the plane fitted to the other polygon.
  bool nearly_coplanar(const std::vector<Point_3>& a, const Plane_3& plane_a, const std::vector<Point_3>& b, const Plane_3& plane_b) const {
    // Relative tolerance. Smaller deviations lead to very thin volumes in the kinetic partition.
    const FT tolerance = FT(1e-6);

    FT min[3], max[3];
    for (std::size_t d = 0; d < 3; d++) {
      min[d] = max[d] = a[0][d];
      for (const Point_3& p : a) {
        min[d] = (std::min)(min[d], p[d]);
        max[d] = (std::max)(max[d], p[d]);
      }
      for (const Point_3& p : b) {
        min[d] = (std::min)(min[d], p[d]);
        max[d] = (std::max)(max[d], p[d]);
      }
    }
    const FT sq_extent = CGAL::square(max[0] - min[0]) + CGAL::square(max[1] - min[1]) + CGAL::square(max[2] - min[2]);
    const FT sq_tolerance = CGAL::square(tolerance) * sq_extent;

    // Squared distance of p to the plane without requiring a distance functor of the traits.
    auto sq_distance = [](const Point_3& p, const Plane_3& pl) {
      const FT v = pl.a() * p.x() + pl.b() * p.y() + pl.c() * p.z() + pl.d();
      return CGAL::square(v) / (CGAL::square(pl.a()) + CGAL::square(pl.b()) + CGAL::square(pl.c()));
    };

    for (const Point_3& p : a)
      if (sq_distance(p, plane_b) > sq_tolerance)
        return false;

    for (const Point_3& p : b)
      if (sq_distance(p, plane_a) > sq_tolerance)
        return false;

    return true;
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
        CGAL_assertion(false);
      }

      if (!pos && !neg) {
        std::cout << "face is degenerated" << std::endl;
        CGAL_assertion(false);
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
        CGAL_assertion(pts_idx[f][v].first != static_cast<std::size_t>(-1));
        CGAL_assertion(pts_idx[f][v].second != static_cast<std::size_t>(-1));
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
    //std::size_t newpts = 0;
    for (typename CDTplus::Finite_vertices_iterator vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit) {
      if (!vit->info().input) {
        vit->info().point_3 = plane.to_3d(vit->point());
        vit->info().idA2 = vit->info().idB2 = Index(-1, -1);
        //newpts++;
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

    //std::size_t newpts = 0;
    // Generate 3D points corresponding to the intersections
    for (typename CDTplus::Finite_vertices_iterator vit = cdtC.finite_vertices_begin(); vit != cdtC.finite_vertices_end(); ++vit) {
      if (!vit->info().input) {
        vit->info().point_3 = plane.to_3d(vit->point());
        vit->info().idA2 = vit->info().idB2 = Index(-1, -1);
        //newpts++;
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

  // Discards the result of a previous partition.
  void clear_partition() {
    m_volumes.clear();
    m_index2volume.clear();
  }

  // Releases the kinetic data structures kept by initialize().
  void release_kinetic_data() {
    for (Sub_partition& partition : m_partition_nodes)
      partition.m_data.reset();

    m_kinetic_data_initialized = false;
  }

  // Sets the parameters and separates the input polygons into the leaves of an octree.
  template<typename NamedParameters>
  bool create_subpartitions(const NamedParameters& np) {
    Timer timer;
    m_parameters.bbox_dilation_ratio = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::bbox_dilation_ratio), FT(11) / FT(10));
    m_parameters.reorient_bbox = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::reorient_bbox), false);
    m_parameters.max_octree_depth = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::max_octree_depth), 3);
    m_parameters.max_octree_node_size = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::max_octree_node_size), 40);

    if (m_parameters.bbox_dilation_ratio <= 1.0) {
      std::cout << "Warning: You set enlarge_bbox_ratio <= 1.0! The valid range is (1.0, +inf). Setting to 1.1!" << std::endl;
      m_parameters.bbox_dilation_ratio = FT(11) / FT(10);
    }

    std::cout.precision(17);
    if (m_input_polygons.size() == 0) {
      std::cout << "Warning: Your input is empty!" << std::endl;
      return false;
    }

    std::set<std::size_t> n;
    for (auto p : m_input2regularized)
      if (p != std::size_t(-1))
        n.insert(p);

    CGAL_assertion(m_regularized2input.size() == m_input_polygons.size());
    CGAL_assertion(m_regularized2input.size() == n.size());

    if (m_parameters.bbox_dilation_ratio < FT(1)) {
      CGAL_warning_msg(m_parameters.bbox_dilation_ratio >= FT(1),
        "Warning: You set enlarge_bbox_ratio < 1.0! The valid range is [1.0, +inf). Setting to 1.0!");
      m_parameters.bbox_dilation_ratio = FT(1);
    }

    if (m_parameters.verbose) {
      //const std::string is_reorient = (m_parameters.reorient ? "true" : "false");

      std::cout << "\n--- PARTITION OPTIONS: \n";
      std::cout << "* enlarge bbox ratio: " << m_parameters.bbox_dilation_ratio << "\n";

      std::cout << "\n--- CREATING SUBPARTITIONS:" << std::endl;

      // Initialization.
      timer.reset();
      timer.start();
    }

    if (m_parameters.debug) {
      for (std::size_t i = 0; i < m_input_polygons.size(); i++)
        KSP_3::internal::dump_polygon<GeomTraits>(m_input_polygons[i], std::to_string(i) + "-input_polygon");
    }

    if (!split_octree())
      return false;
    m_partitions.resize(m_partition_nodes.size());
    std::iota(m_partitions.begin(), m_partitions.end(), 0);

    for (std::size_t idx : m_partitions)
      m_partition_nodes[idx].index = idx;

    // Timing.
    if (m_parameters.verbose) {
      const double time_to_initialize = timer.time();
      std::cout << "* octree creation time: " << time_to_initialize << std::endl;
    }

    return true;
  }

  // Calls f(i) for each octree leaf m_partitions[i], in parallel if ConcurrencyTag is Parallel_tag.
  // The leaves are processed in the order of decreasing number of input polygons, so that the most expensive leaves are
  // started first and do not delay the end of the computation. Each leaf is scheduled separately for the same reason.
  template<typename ConcurrencyTag, typename Function>
  void for_each_leaf(const Function& f) const {
    std::vector<std::size_t> order(m_partitions.size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
      return m_partition_nodes[m_partitions[a]].input_polygons.size() > m_partition_nodes[m_partitions[b]].input_polygons.size();
    });

#ifdef CGAL_LINKED_WITH_TBB
    if constexpr (std::is_convertible_v<ConcurrencyTag, Parallel_tag>) {
      tbb::parallel_for(tbb::blocked_range<std::size_t>(0, order.size(), 1), [&](const tbb::blocked_range<std::size_t>& r) {
        for (std::size_t i = r.begin(); i != r.end(); ++i)
          f(order[i]);
      }, tbb::simple_partitioner());
      return;
    }
#endif

    for (std::size_t i : order)
      f(i);
  }

  // Initializes the kinetic data structures of all octree leaves.
  template<typename ConcurrencyTag>
  bool initialize_subpartitions(FT& initialization_time) {
    // Times are collected per partition. In the parallel version, they sum up to more than the elapsed time.
    std::vector<FT> initialization_times(m_partitions.size(), 0);
    std::atomic<bool> failed(false);

    for_each_leaf<ConcurrencyTag>([&](const std::size_t i) {
      if (failed)
        return;

      if (!initialize_subpartition(m_partitions[i], initialization_times[i]))
        failed = true;
    });

    initialization_time = std::accumulate(initialization_times.begin(), initialization_times.end(), FT(0));

    return !failed;
  }

  // Computes the kinetic partitions of all octree leaves.
  // If init_and_release is true, the kinetic data structure of each leaf is created before and released after its propagation.
  // Otherwise, the kinetic data structures created by initialize() are used and restored afterwards.
  template<typename ConcurrencyTag>
  bool compute_subpartitions(std::size_t k, bool init_and_release, FT& initialization_time, FT& partition_time, FT& finalization_time) {
    // Times are collected per partition. In the parallel version, they sum up to more than the elapsed time.
    std::vector<FT> initialization_times(m_partitions.size(), 0), partition_times(m_partitions.size(), 0), finalization_times(m_partitions.size(), 0);
    std::atomic<bool> failed(false);

    std::cout.precision(17);
    if (m_parameters.verbose) {
      if (std::is_same<ConcurrencyTag, Sequential_tag>::value)
        std::cout << "sequential" << std::endl;
      else if (std::is_same<ConcurrencyTag, Parallel_tag>::value)
        std::cout << "parallel" << std::endl;
    }

    for_each_leaf<ConcurrencyTag>([&](const std::size_t i) {
      if (failed)
        return;

      if (!compute_partition(m_partitions[i], k, init_and_release, initialization_times[i], partition_times[i], finalization_times[i])) {
        failed = true;
        KSP::internal::Log_line() << "Warning: partition() failed for partition " << i << " " << m_partitions[i];
      }
    });

    initialization_time = std::accumulate(initialization_times.begin(), initialization_times.end(), FT(0));
    partition_time = std::accumulate(partition_times.begin(), partition_times.end(), FT(0));
    finalization_time = std::accumulate(finalization_times.begin(), finalization_times.end(), FT(0));

    return !failed;
  }

  // Creates and initializes the kinetic data structure of one octree leaf.
  // Only accesses data of the given Sub_partition and read-only shared data. Thus, it can be called in parallel for different leaves.
  bool initialize_subpartition(std::size_t idx, FT& initialization_time) {
    Sub_partition& partition = m_partition_nodes[idx];
    // Each line written to std::cout gets the index of the partition as prefix, as partitions may be computed in parallel.
    const std::string log_prefix = "[P" + std::to_string(idx) + "] ";
    Timer timer;
    timer.start();

    partition.m_data = std::make_shared<Data_structure>(std::to_string(idx) + "-", log_prefix);

    Initializer initializer(partition.clipped_polygons, partition.m_input_planes, *partition.m_data, m_parameters, m_partitions.size() != 1);
    initializer.initialize(partition.bbox, partition.input_polygons);

    // The input of the Initializer is not needed anymore.
    std::vector<std::vector<typename Intersection_kernel::Point_3> >().swap(partition.clipped_polygons);
    std::vector<typename Intersection_kernel::Plane_3>().swap(partition.m_input_planes);

    initialization_time = timer.time();

    // Initialization successful?
    if (partition.m_data->number_of_support_planes() < 6) {
      KSP::internal::Log_line(log_prefix) << "Kinetic partition not initialized or empty. Number of support planes: " << partition.m_data->number_of_support_planes();

      return false;
    }

    return true;
  }

  // Propagates and finalizes the kinetic partition of one octree leaf. If init_and_release is true, the kinetic data structure is
  // created before and released afterwards. Otherwise, the kinetic data structure is restored to its initialized state afterwards.
  // Only accesses data of the given Sub_partition and read-only shared data. Thus, it can be called in parallel for different leaves.
  bool compute_partition(std::size_t idx, std::size_t k, bool init_and_release, FT& initialization_time, FT& partition_time, FT& finalization_time) {
    Sub_partition& partition = m_partition_nodes[idx];
    const std::string log_prefix = "[P" + std::to_string(idx) + "] ";

    if (init_and_release) {
      if (!initialize_subpartition(idx, initialization_time))
        return false;
    }

    CGAL_assertion(partition.m_data != nullptr);

    // Propagation and finalization modify the kinetic data structure. Its initialized state is saved to be restored afterwards.
    typename Data_structure::Initial_state initial_state;
    if (!init_and_release)
      initial_state = partition.m_data->save_initial_state();

    Timer timer;
    timer.start();

    if (m_parameters.verbose) {
      KSP::internal::Log_line(log_prefix) << "--- RUNNING THE QUEUE:";
      KSP::internal::Log_line(log_prefix) << "* propagation started";
    }

    // Propagation.
    Propagation propagation(*partition.m_data, m_parameters);
    propagation.propagate(k);

    partition_time = timer.time();
    timer.reset();

    if (m_parameters.verbose) {
      KSP::internal::Log_line(log_prefix) << "events in partition: " << propagation.handled_events;
      KSP::internal::Log_line(log_prefix) << "* propagation finished";
      KSP::internal::Log_line(log_prefix) << "--- FINALIZING PARTITION:";
    }

    // Finalization.
    for (std::size_t i = 0; i < partition.m_data->number_of_support_planes(); i++)
      if (!partition.m_data->support_plane(i).mesh().is_valid(true))
        KSP::internal::Log_line(log_prefix) << i << ". support has an invalid mesh!";

    Finalizer finalizer(*partition.m_data, m_parameters);

    if (m_parameters.verbose)
      KSP::internal::Log_line(log_prefix) << "* getting volumes ...";

    finalizer.create_polyhedra();
    finalization_time = timer.time();

    if (m_parameters.verbose)
      KSP::internal::Log_line(log_prefix) << "partition with " << partition.input_polygons.size() << " input polygons split into " << partition.m_data->number_of_volumes() << " volumes";

    extract_finalized_data(partition, idx);

    if (init_and_release) {
      // The kinetic data structure (support planes, intersection graph, kinetic vertices, ...) is not needed anymore.
      // Releasing it here keeps the peak memory at one kinetic data structure per thread plus the finalized results.
      CGAL_assertion(partition.m_data.use_count() == 1);
      partition.m_data.reset();
    }
    else
      partition.m_data->restore_initial_state(std::move(initial_state));

    return true;
  }

  // Moves the results of the Finalizer from the kinetic data structure into the Sub_partition.
  // Afterwards, make_conformal() and all access functions only use the data in the Sub_partition.
  void extract_finalized_data(Sub_partition& partition, std::size_t idx) {
    Data_structure& data = *partition.m_data;

    // Volumes
    partition.volumes.resize(data.volumes().size());
    for (std::size_t i = 0; i < data.volumes().size(); i++) {
      partition.volumes[i].faces = std::move(data.volumes()[i].faces);
      partition.volumes[i].centroid = data.volumes()[i].centroid;
    }

    // Vertices
    partition.vertices = std::move(data.vertices());
    partition.exact_vertices = std::move(data.exact_vertices());

    // Faces
    partition.face_neighbors.resize(data.face_to_volumes().size());
    for (std::size_t j = 0; j < data.face_to_volumes().size(); j++) {
      const auto& p = data.face_to_volumes()[j];
      partition.face_neighbors[j] = std::make_pair(Index(idx, p.first), Index(idx, p.second));
    }

    partition.face2vertices.resize(data.face_to_vertices().size());
    for (std::size_t j = 0; j < data.face_to_vertices().size(); j++) {
      const auto& f = data.face_to_vertices()[j];
      partition.face2vertices[j].resize(f.size());
      for (std::size_t k = 0; k < f.size(); k++)
        partition.face2vertices[j][k] = std::make_pair(idx, f[k]);
    }

    partition.face2sp = std::move(data.face_to_support_plane());
    partition.face_part_of_input_polygon = std::move(data.face_is_part_of_input_polygon());
    partition.face_input = std::move(data.face_to_input_polygon());

    // Support planes
    partition.sp_exact_planes.resize(data.number_of_support_planes());
    partition.sp_input_polygon.resize(data.number_of_support_planes());
    for (std::size_t i = 0; i < data.number_of_support_planes(); i++) {
      partition.sp_exact_planes[i] = data.support_plane(i).exact_plane();
      partition.sp_input_polygon[i] = data.support_plane(i).data().actual_input_polygon;
    }

    partition.input2sp = data.input_polygon_map();
  }

  // Merges the finalized partitions of the octree leaves into one conformal partition.
  void merge_subpartitions(FT& conformal_time) {
    for (std::size_t i = 0; i < m_partitions.size(); i++) {
      const Sub_partition& partition = m_partition_nodes[m_partitions[i]];

      for (std::size_t j = 0; j < partition.volumes.size(); j++)
        m_volumes.push_back(std::make_pair(m_partitions[i], j));
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
          auto it = pts2idx.emplace(m_partition_nodes[face[j].first].exact_vertices[face[j].second], face[j]);
          if (!it.second)
            face[j] = it.first->second;
        }
      }
    }

    Timer timer;
    timer.start();
    make_conformal(0);
    conformal_time = timer.time();

    if (m_parameters.verbose)
      std::cout << "ksp v: " << m_partition_nodes[0].vertices.size() << " f: " << m_partition_nodes[0].face2vertices.size() << " vol: " << m_volumes.size() << std::endl;
  }

  void collect_faces(std::size_t partition_idx, std::size_t sp_idx, std::vector<Index>& faces, typename Intersection_kernel::Plane_3& plane) {
    Sub_partition& p = m_partition_nodes[partition_idx];

    plane = p.sp_exact_planes[sp_idx];

    const std::vector<std::size_t>& f2sp = p.face2sp;

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

  template<typename VertexMap>
  bool can_add_face_to_lcc(const std::vector<Index>& vtx_of_face, const VertexMap& mapped_vertices,
                           std::vector<bool>& used_vertices) const {
    auto it = mapped_vertices.find(vtx_of_face.back());
    CGAL_assertion(it != mapped_vertices.end());
    bool before = used_vertices[it->second];

    for (std::size_t i = 0; i < vtx_of_face.size(); i++) {
      it = mapped_vertices.find(vtx_of_face[i]);
      CGAL_assertion(it != mapped_vertices.end());
      if (used_vertices[it->second]) {
        if (before)
          return true;
        else
          before = true;
      }
      else
        before = false;
    }

    return false;
  }

  template<typename LCC, typename LCCIncrementalBuilder, typename VertexMap>
  void add_face_to_lcc(LCC &lcc, LCCIncrementalBuilder &ib, const Index& face, std::vector<Index> &vtx_of_face,
                       const std::vector<std::size_t>& vtx_of_volume, VertexMap &mapped_vertices, const std::vector<IkPoint_3> &vtx,
                       std::vector<bool> &used_vertices) const {
    ib.begin_facet();

    // The orientation of the face is decided exactly, which is required for thin volumes, where an approximated centroid may lie outside.
    // The points are chosen with their approximations, such that the exact predicates are usually decided by their filters.
    auto approx = [&](std::size_t vi) {
      return std::array<double, 3>{ CGAL::to_double(vtx[vi].x()), CGAL::to_double(vtx[vi].y()), CGAL::to_double(vtx[vi].z()) };
    };
    auto cross = [](const std::array<double, 3>& a, const std::array<double, 3>& b, const std::array<double, 3>& c) {
      const double u[3] = { b[0] - a[0], b[1] - a[1], b[2] - a[2] }, v[3] = { c[0] - a[0], c[1] - a[1], c[2] - a[2] };
      return std::array<double, 3>{ u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0] };
    };

    // Three consecutive vertices of the face spanning the largest approximated area. The face is convex, so they are not collinear.
    const std::size_t nv = vtx_of_face.size();
    std::size_t best = 0;
    double best_area = -1;
    for (std::size_t k = 0; k < nv; k++) {
      const auto c = cross(approx(mapped_vertices[vtx_of_face[k]]), approx(mapped_vertices[vtx_of_face[(k + 1) % nv]]), approx(mapped_vertices[vtx_of_face[(k + 2) % nv]]));
      const double area = c[0] * c[0] + c[1] * c[1] + c[2] * c[2];
      if (area > best_area) {
        best_area = area;
        best = k;
      }
    }
    if (CGAL::collinear(vtx[mapped_vertices[vtx_of_face[best]]], vtx[mapped_vertices[vtx_of_face[(best + 1) % nv]]], vtx[mapped_vertices[vtx_of_face[(best + 2) % nv]]])) {
      // The approximation failed, search exactly.
      for (best = 0; best < nv; best++)
        if (!CGAL::collinear(vtx[mapped_vertices[vtx_of_face[best]]], vtx[mapped_vertices[vtx_of_face[(best + 1) % nv]]], vtx[mapped_vertices[vtx_of_face[(best + 2) % nv]]]))
          break;
    }
    CGAL_assertion(best < nv);
    const std::size_t pi = mapped_vertices[vtx_of_face[best]], qi = mapped_vertices[vtx_of_face[(best + 1) % nv]], ri = mapped_vertices[vtx_of_face[(best + 2) % nv]];
    const IkPoint_3& p = vtx[pi];
    const IkPoint_3& q = vtx[qi];
    const IkPoint_3& r = vtx[ri];

    // The volume is convex. Thus, every vertex of the volume that is not on the supporting plane of the face lies on its inner side.
    // The vertex with the largest approximated distance to the plane is tested first.
    const std::array<double, 3> ap = approx(pi), normal = cross(ap, approx(qi), approx(ri));
    std::size_t far_vi = vtx_of_volume.front();
    double far_dist = -1;
    for (std::size_t vi : vtx_of_volume) {
      const std::array<double, 3> a = approx(vi);
      const double d = std::abs(normal[0] * (a[0] - ap[0]) + normal[1] * (a[1] - ap[1]) + normal[2] * (a[2] - ap[2]));
      if (d > far_dist) {
        far_dist = d;
        far_vi = vi;
      }
    }
    CGAL::Orientation inner = CGAL::orientation(p, q, r, vtx[far_vi]);
    if (inner == CGAL::COPLANAR) {
      // The approximation failed, search exactly.
      for (std::size_t vi : vtx_of_volume) {
        inner = CGAL::orientation(p, q, r, vtx[vi]);
        if (inner != CGAL::COPLANAR)
          break;
      }
    }
    CGAL_assertion(inner != CGAL::COPLANAR);

    // Faces are inserted such that the volume lies on the positive side of their vertex order.
    if (inner != CGAL::POSITIVE)
      std::reverse(vtx_of_face.begin(), vtx_of_face.end());

    for (const Index& v : vtx_of_face) {
      ib.add_vertex_to_facet(static_cast<typename LCC::size_type>(mapped_vertices[v]));
      if (!used_vertices[mapped_vertices[v]]) {
        used_vertices[mapped_vertices[v]] = true;
      }
    }

    auto face_dart = ib.end_facet(); // returns a dart to the face
    if (lcc.template attribute<2>(face_dart) == lcc.null_descriptor) {
      lcc.template set_attribute<2>(face_dart, lcc.template create_attribute<2>());
      // How to handle bbox planes that coincide with input polygons? Check support plane
      std::size_t sp = m_partition_nodes[face.first].face2sp[face.second];

      // There are three different cases:
      // 1. face belongs to a plane from an input polygon
      // 2. face originates from octree splitting (and does not have an input plane)
      // 3. face lies on the bbox
      // Faces on a support plane of several coplanar input polygons refer to the input polygon they originate from.
      // Faces on the bbox may coincide with an input polygon of the support plane.
      int ip = static_cast<int>(m_partition_nodes[face.first].face_input[face.second]);
      if (ip == -1)
        ip = static_cast<int>(m_partition_nodes[face.first].sp_input_polygon[sp]);

      if (ip != -1) {
        // Index of the inserted polygon. It differs from the index of the input polygon of the partition if degenerate polygons have been dropped.
        const std::size_t input = m_regularized2input[m_partition_nodes[face.first].input_polygons[ip]].front();
        lcc.template info<2>(face_dart).input_polygon_index = static_cast<Face_support>(input);
        lcc.template info<2>(face_dart).input_plane_index = static_cast<int>(m_input_plane_indices[input]);
      }
      else {
        // If there is no input polygon, check whether it has two neighbors
        auto n = neighbors(face);
        if (n.second >= 0)
          lcc.template info<2>(face_dart).input_polygon_index = static_cast<Face_support>(-7);
        else
          lcc.template info<2>(face_dart).input_polygon_index = static_cast<Face_support>(n.second);
        lcc.template info<2>(face_dart).input_plane_index = lcc.template info<2>(face_dart).input_polygon_index;
      }

      lcc.template info<2>(face_dart).part_of_initial_polygon = m_partition_nodes[face.first].face_part_of_input_polygon[face.second];

      lcc.template info<2>(face_dart).plane = m_partition_nodes[face.first].sp_exact_planes[sp];
    }
    else if (!lcc.template info<2>(face_dart).part_of_initial_polygon && m_partition_nodes[face.first].face_part_of_input_polygon[face.second]) {
      // A face on a split plane of the octree is present in both adjacent leaves. An input polygon lying in the split plane
      // may only be assigned to one of the leaves, so only one of them marks the face as part of the input polygon.
      const int ip = static_cast<int>(m_partition_nodes[face.first].face_input[face.second]);
      if (ip != -1) {
        const std::size_t input = m_regularized2input[m_partition_nodes[face.first].input_polygons[ip]].front();
        lcc.template info<2>(face_dart).input_polygon_index = static_cast<Face_support>(input);
        lcc.template info<2>(face_dart).input_plane_index = static_cast<int>(m_input_plane_indices[input]);
      }
      lcc.template info<2>(face_dart).part_of_initial_polygon = true;
    }
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
    FT la = CGAL::approximate_sqrt(axis1.squared_length());
    axis1 = axis1 * (1.0 / la);
    FT lb = CGAL::approximate_sqrt(axis2.squared_length());
    axis2 = axis2 * (1.0 / lb);

    if (CGAL::abs(axis1.x()) < CGAL::abs(axis2.x())) {
      Vector_2 tmp = axis1;
      axis1 = axis2;
      axis2 = tmp;
    }

    if (0 > axis1.x())
      axis1 = -axis1;

    axis2 = Vector_2(-axis1.y(), axis1.x());

    CGAL::Aff_transformation_3<Kernel> R(axis1.x(), axis1.y(), 0,
                                        -axis1.y(), axis1.x(), 0,
                                                 0,         0, 1.0);

    CGAL::Aff_transformation_3<Kernel> T(CGAL::TRANSLATION, -Vector_3((bbox[0].x() + bbox[2].x()) * 0.5, (bbox[0].y() + bbox[2].y()) * 0.5, (maxz + minz) * 0.5));

    return R * T;
  }

  bool same_face(const Face_handle& a, const Face_handle& b) const {
    return (b->info().idA2 == a->info().idA2 && b->info().idB2 == a->info().idB2);
  }

  void set_face(const Index& f, const Index& other, std::set<Index>& replaced, const std::vector<Vertex_handle>& polygon) {
    From_exact from_exact;
    auto pair = replaced.insert(f);
    std::size_t idx;
    CGAL_assertion(m_partition_nodes[f.first].face_neighbors[f.second].first.first == f.first);
    std::size_t vol_idx = m_partition_nodes[f.first].face_neighbors[f.second].first.second;

    if (!pair.second) {
      // New face has a new index
      idx = m_partition_nodes[f.first].face2vertices.size();
      // Add face into vector
      m_partition_nodes[f.first].face2vertices.push_back(std::vector<Index>());
      m_partition_nodes[f.first].face_part_of_input_polygon.push_back(m_partition_nodes[f.first].face_part_of_input_polygon[f.second]);
      // Add face index into volume
      m_partition_nodes[f.first].volumes[vol_idx].faces.push_back(idx);
      // Copy neighbor from already existing face
      m_partition_nodes[f.first].face_neighbors.push_back(m_partition_nodes[f.first].face_neighbors[f.second]);
      m_partition_nodes[f.first].face2sp.push_back(m_partition_nodes[f.first].face2sp[f.second]);
      m_partition_nodes[f.first].face_input.push_back(m_partition_nodes[f.first].face_input[f.second]);
    }
    else {
      idx = f.second;
      // All boundary faces should have a negative second neighbor.
      CGAL_assertion(m_partition_nodes[f.first].face_neighbors[idx].second.second >= std::size_t(-6));
    }
    std::vector<Index>& vertices = m_partition_nodes[f.first].face2vertices[idx];
    // First neighbor of other face should point to the inside volume in the other partition and thus cannot be negative
    CGAL_assertion(m_partition_nodes[other.first].face_neighbors[other.second].first.second < std::size_t(-6));
    m_partition_nodes[f.first].face_neighbors[idx].second = m_partition_nodes[other.first].face_neighbors[other.second].first;
    vertices.resize(polygon.size());
    for (std::size_t i = 0; i < polygon.size(); i++) {
      VI& vi = polygon[i]->info();

      if (vi.idA2.first < vi.idB2.first)
        vertices[i] = vi.idA2;
      else if (vi.idB2.first != static_cast<std::size_t>(-1))
        vertices[i] = vi.idB2;
      else {
        std::size_t vidx = m_partition_nodes[f.first].vertices.size();
        m_partition_nodes[f.first].vertices.push_back(from_exact(vi.point_3));
        m_partition_nodes[f.first].exact_vertices.push_back(vi.point_3);
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
        Face_handle last = fh;

        // Mark seed face as segmented
        fh->info().id2.first = extracted;

        // edge is pair<Face_handle, int (vertex)>
        while (face.front() != face.back()) {
          auto eit = cdt.incident_edges(face.back(), last);
          auto first = eit;

          CGAL_assertion(!cdt.is_infinite(eit->first));
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
            CGAL_assertion(eit != first);
          } while (eit != first);

          CGAL_assertion(face.size() < 100);
        }

        // The last vertex is equal to the first one, so it should be removed.
        face.pop_back();

        // face ids in partitions are stored in fh->info
        ID& id = fh->info();
        set_face(id.idA2, id.idB2, replacedA, face);
        set_face(id.idB2, id.idA2, replacedB, face);
      }
    }
  }

  std::pair<std::size_t, int> find_portal(std::size_t volume, std::size_t former, const Index& vA, const Index& vB, std::size_t& portal) const {
    portal = static_cast<std::size_t>(-7);
    auto vol = m_volumes[volume];
    const std::vector<std::size_t>& faces = m_partition_nodes[vol.first].volumes[vol.second].faces;

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

  void adapt_internal_edges(const CDTplus& cdtC, std::vector<std::vector<Constraint_info> >& c) {
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

        // Extract edge
        std::vector<Index> vertices_of_edge;
        for (typename CDTplus::Vertices_in_constraint_iterator vi = cdtC.vertices_in_constraint_begin(id); vi != cdtC.vertices_in_constraint_end(id); vi++) {
          if ((*vi)->info().idA2.first == static_cast<std::size_t>(-1))
            vertices_of_edge.push_back((*vi)->info().idB2);
          else vertices_of_edge.push_back((*vi)->info().idA2);
        }

        if (vertices_of_edge.size() == 2)
          continue;

        // Check length of constraint
        // size 2 means it has not been split, thus there are no t-junctions.
        CGAL_assertion(vertices_of_edge.size() >= 2);

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
    std::map<std::size_t, std::vector<Index> > a_sets, b_sets;
    for (const Index& i : a)
      a_sets[i.first].push_back(i);
    for (const Index& i : b)
      b_sets[i.first].push_back(i);

    // At first, one CDTplus is created for each partition node
    std::vector<CDTplus> a_cdts(a_sets.size()), b_cdts(b_sets.size());

    std::vector< std::vector<std::vector<Constraint_info> > > a_constraints;
    std::vector< std::vector<std::vector<Constraint_info> > > b_constraints;


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

    for (idx = 0; idx < a_sets.size(); idx++)
      adapt_internal_edges(cdtC, a_constraints[idx]);

    for (idx = 0; idx < b_sets.size(); idx++)
      adapt_internal_edges(cdtC, b_constraints[idx]);
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
    for (std::size_t dim = 0; dim < 3; dim++) {
      std::vector<Index> lower, upper;
      typename Intersection_kernel::Plane_3 plane;

      collect_opposing_faces(node, dim, lower, upper, plane);

      make_conformal(lower, upper, plane);
    }
  }

  // Projects p onto the plane along the coordinate axis that is most orthogonal to the plane, i.e., only the coordinate
  // of that axis is changed. Like an orthogonal projection, the point is afterwards exactly on the plane, but the other
  // two coordinates are kept as they are. This keeps the bit complexity of the exact coordinates low.
  static typename Intersection_kernel::Point_3 project_along_axis(const typename Intersection_kernel::Plane_3& pl, const typename Intersection_kernel::Point_3& p) {
    using IFT = typename Intersection_kernel::FT;
    const IFT& a = pl.a();
    const IFT& b = pl.b();
    const IFT& c = pl.c();
    const IFT& d = pl.d();
    const double ax = CGAL::abs(CGAL::to_double(a)), ay = CGAL::abs(CGAL::to_double(b)), az = CGAL::abs(CGAL::to_double(c));
    // The approximation of the largest coefficient is non-zero, thus the coefficient itself is non-zero.
    if (az >= ax && az >= ay)
      return typename Intersection_kernel::Point_3(p.x(), p.y(), -(a * p.x() + b * p.y() + d) / c);
    if (ay >= ax)
      return typename Intersection_kernel::Point_3(p.x(), -(a * p.x() + c * p.z() + d) / b, p.z());
    return typename Intersection_kernel::Point_3(-(b * p.y() + c * p.z() + d) / a, p.y(), p.z());
  }

  // Returns false if the input polygons do not span a volume, e.g., if all of them lie in the same plane.
  bool split_octree() {
    // Octree creation for sub partition
    std::size_t count = 0;
    for (const auto& p : m_input_polygons)
      count += p.size();

    m_points.clear();
    m_points.reserve(count);
    m_polygons.clear();
    m_polygons.reserve(m_input_polygons.size());

    To_exact to_exact;
    From_exact from_exact;

    std::vector<typename Intersection_kernel::Point_3> pts;
    pts.reserve(count);

    if (m_parameters.reorient_bbox) {
      m_transform = to_exact(get_obb2abb(m_input_polygons));
      auto trans = from_exact(m_transform);

      for (std::size_t i = 0; i < m_input_polygons.size(); i++) {
        const auto& p = m_input_polygons[i];
        std::size_t idx = m_points.size();
        for (const Point_3& pt : p) {
          m_points.push_back(trans.transform(pt));
          // Project exactly onto the support plane, so that the clipped polygons are exactly
          // incident to the intersection lines of their support plane and the octree box planes.
          pts.push_back(m_transform.transform(project_along_axis(m_input_planes[i], to_exact(pt))));
        }

        m_polygons.push_back(std::vector<std::size_t>(p.size()));
        std::iota(m_polygons.back().begin(), m_polygons.back().end(), idx);
      }
    }
    else {
      m_transform = CGAL::Aff_transformation_3<Intersection_kernel>(CGAL::IDENTITY);

      for (std::size_t i = 0; i < m_input_polygons.size(); i++) {
        const auto& p = m_input_polygons[i];
        std::size_t idx = m_points.size();
        std::copy(p.begin(), p.end(), std::back_inserter(m_points));
        // Project exactly onto the support plane, so that the clipped polygons are exactly
        // incident to the intersection lines of their support plane and the octree box planes.
        for (const auto &pt : p)
          pts.push_back(project_along_axis(m_input_planes[i], to_exact(pt)));
        m_polygons.push_back(std::vector<std::size_t>(p.size()));
        std::iota(m_polygons.back().begin(), m_polygons.back().end(), idx);
      }
    }

    // The bounding box, which is dilated per dimension, needs a positive extent in each dimension.
    for (std::size_t d = 0; d < 3; d++) {
      auto cmp = [d](const typename Intersection_kernel::Point_3& a, const typename Intersection_kernel::Point_3& b) { return a[d] < b[d]; };
      auto mm = std::minmax_element(pts.begin(), pts.end(), cmp);
      if (mm.first == pts.end() || (*mm.first)[d] == (*mm.second)[d]) {
        KSP::internal::Log_line() << "Warning: The input polygons lie in a plane " << (m_parameters.reorient_bbox ? "" : "orthogonal to an axis ")
          << "and do not span a volume. No partition is created.";
        return false;
      }
    }

    m_octree = std::make_unique<Octree>(CGAL::Orthtree_traits_polygons<Intersection_kernel>(pts, m_polygons, m_parameters.bbox_dilation_ratio));
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

    m_partition_nodes.clear();
    m_partition_nodes.resize(leaf_count);

    m_node2partition.resize(max_count + 1, std::size_t(-1));

    std::size_t idx = 0;
    CGAL::Aff_transformation_3<Intersection_kernel> inv = m_transform.inverse();
    for (typename Octree::Node_index node : m_octree->traverse(CGAL::Orthtrees::Leaves_traversal<Octree>(*m_octree)))
      if (m_octree->is_leaf(node)) {
        // Creating bounding box
        CGAL::Iso_cuboid_3<Intersection_kernel> box = m_octree->bbox(node);
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

        // The octree clipping assigns a polygon that only touches the box of the node in a vertex or along an edge
        // to that node as well.
        polys.erase(std::remove_if(polys.begin(), polys.end(), [](const auto& p) {
          const auto& pts = p.second;
          std::size_t j = 1;
          while (j < pts.size() && pts[j] == pts[0])
            j++;
          for (std::size_t k = j + 1; k < pts.size(); k++)
            if (!CGAL::collinear(pts[0], pts[j], pts[k]))
              return false;
          return true;
          }), polys.end());

        for (std::size_t j = 0; j < polys.size(); j++) {
          m_partition_nodes[idx].input_polygons.push_back(polys[j].first);
          m_partition_nodes[idx].m_input_planes.push_back(m_input_planes[polys[j].first]);
        }

        m_partition_nodes[idx].clipped_polygons.resize(polys.size());
        //auto inv_non_exact = from_exact(inv);
        for (std::size_t i = 0; i < polys.size(); i++) {
          m_partition_nodes[idx].clipped_polygons[i].resize(polys[i].second.size());
          for (std::size_t j = 0; j < polys[i].second.size(); j++) {
            m_partition_nodes[idx].clipped_polygons[i][j] = inv.transform(polys[i].second[j]);
          }
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

          //KSP_3::internal::dump_polygons(m_partition_nodes[idx].clipped_polygons, std::to_string(idx) + "-polys.ply");
        }
        idx++;
      }
    if (m_parameters.verbose)
      std::cout << "input split into " << m_partition_nodes.size() << " partitions" << std::endl;

    return true;
  }
};

inline
void write_cmap_attribute_node(boost::property_tree::ptree& node, const typename Kinetic_space_partition_3<CGAL::Epick>::Linear_cell_complex_min_items::Face_attribute& fa) {
  node.add("v.a", fa.plane.a());
  node.add("v.b", fa.plane.b());
  node.add("v.c", fa.plane.c());
  node.add("v.d", fa.plane.d());
  node.add("v.ip", fa.input_polygon_index);
  node.add("v.ipl", fa.input_plane_index);
  node.add("v.poip", fa.part_of_initial_polygon);
}

inline
void write_cmap_attribute_node(boost::property_tree::ptree& node, const typename Kinetic_space_partition_3<CGAL::Epick>::Linear_cell_complex_min_items::Volume_attribute& va) {
  node.add("v.x", va.barycenter.x());
  node.add("v.y", va.barycenter.y());
  node.add("v.z", va.barycenter.z());
  node.add("v.id", va.volume_id);
}

template<> inline
void read_cmap_attribute_node
(const boost::property_tree::ptree::value_type& v, Exact_predicates_exact_constructions_kernel::Point_3& val) {
  double x = v.second.get<double>("x");
  double y = v.second.get<double>("y");
  double z = v.second.get<double>("z");
  val = Exact_predicates_exact_constructions_kernel::Point_3(x, y, z);
}

template<> inline
void read_cmap_attribute_node
(const boost::property_tree::ptree::value_type& v, typename Kinetic_space_partition_3<Epick>::Linear_cell_complex_min_items::Face_attribute& fa) {
  double a = v.second.get<double>("a");
  double b = v.second.get<double>("b");
  double c = v.second.get<double>("c");
  double d = v.second.get<double>("d");
  fa.plane = Exact_predicates_exact_constructions_kernel::Plane_3(a, b, c, d);
  fa.input_polygon_index = v.second.get<int>("ip");
  fa.input_plane_index = v.second.get<int>("ipl", fa.input_polygon_index);
  fa.part_of_initial_polygon = v.second.get<bool>("poip");
}

template<> inline
void read_cmap_attribute_node
(const boost::property_tree::ptree::value_type& v, typename Kinetic_space_partition_3<Epick>::Linear_cell_complex_min_items::Volume_attribute& va) {
  double x = v.second.get<double>("x");
  double y = v.second.get<double>("y");
  double z = v.second.get<double>("z");
  va.barycenter = Exact_predicates_exact_constructions_kernel::Point_3(x, y, z);
  va.volume_id = v.second.get<int>("id");
}

} // namespace CGAL

#endif // CGAL_KINETIC_SPACE_PARTITION_3_H
