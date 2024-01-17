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

#ifndef CGAL_KINETIC_SURFACE_RECONSTRUCTION_3_H
#define CGAL_KINETIC_SURFACE_RECONSTRUCTION_3_H

#include <CGAL/license/Kinetic_surface_reconstruction.h>

#include <CGAL/Kinetic_space_partition_3.h>
#include <CGAL/KSR_3/Graphcut.h>

#include <CGAL/IO/PLY.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Point_set.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/KSP/debug.h>
#include <CGAL/Shape_regularization/regularize_planes.h>
#include <CGAL/bounding_box.h>

#include <boost/range/adaptor/transformed.hpp>

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/boost/graph/helpers.h>

namespace CGAL
{
/*!
* \ingroup PkgKineticSurfaceReconstructionRef
  \brief Pipeline for piecewise planar surface reconstruction from a point cloud via inside/outside labeling of a kinetic partition using min-cut.

  \tparam GeomTraits
    must be a model of `KineticShapePartitionTraits_3`.

  \tparam PointSet
    must be a range of 3D points and corresponding 3D normal vectors whose iterator type is `RandomAccessIterator`.

  \tparam PointMap
    a model of `ReadablePropertyMap` whose key type is the value type of the `PointSet` and value type is `GeomTraits::Point_3`

  \tparam NormalMap
    a model of `ReadablePropertyMap` whose key type is the value type of the `PointSet` and value type is `GeomTraits::Vector_3`

  \tparam IntersectionKernel
    must be a model of `Kernel` using exact computations. Defaults to `CGAL::Exact_predicates_exact_constructions_kernel`. Used for the internal kinetic shape partition.
*/
template<typename GeomTraits, typename PointRange, typename PointMap, typename NormalMap, typename IntersectionKernel = CGAL::Exact_predicates_exact_constructions_kernel>
class Kinetic_surface_reconstruction_3 {
public:
  using Kernel = GeomTraits;
  using Intersection_kernel = IntersectionKernel;

  using FT = typename Kernel::FT;
  using Point_3 = typename Kernel::Point_3;
  using Plane_3 = typename Kernel::Plane_3;
  using Point_range = PointRange;

  using KSP = Kinetic_space_partition_3<Kernel, Intersection_kernel>;

  using Point_map = PointMap;
  using Normal_map = NormalMap;

  /*!
  \brief Creates a Kinetic_shape_reconstruction_3 object.

  \param points
   an instance of `PointSet` with 3D points and corresponding 3D normal vectors.

  \param np
   a sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

  \cgalNamedParamsBegin
    \cgalParamNBegin{point_map}
      \cgalParamDescription{a property map associating points to the elements of the point set `points`}
      \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type of the iterator of `PointSet` and whose value type is `GeomTraits::Point_3`}
          \cgalParamDefault{`PointMap()`}
    \cgalParamNEnd
  \cgalNamedParamsEnd
  */
  template<typename NamedParameters = parameters::Default_named_parameters>
  Kinetic_surface_reconstruction_3(Point_range& points,
    const NamedParameters& np = CGAL::parameters::default_values()) : m_points(points), m_ground_polygon_index(-1), m_kinetic_partition(np) {
    m_verbose = parameters::choose_parameter(parameters::get_parameter(np, internal_np::verbose), false);
    m_debug = parameters::choose_parameter(parameters::get_parameter(np, internal_np::debug), false);
  }

  /*!
    \brief Detects shapes in the provided point cloud and regularizes them.

    \tparam NamedParameters
    a sequence of \ref bgl_namedparameters "Named Parameters"

    \param np
    an instance of `NamedParameters`.

  \cgalNamedParamsBegin
    \cgalParamNBegin{point_map}
      \cgalParamDescription{a property map associating points to the elements of the point set `points`}
      \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type of the iterator of `PointSet` and whose value type is `GeomTraits::Point_3`}
          \cgalParamDefault{`PointMap()`}
    \cgalParamNEnd
    \cgalParamNBegin{normal_map}
      \cgalParamDescription{a property map associating normals to the elements of the point set `points`}
      \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type of the iterator of `PointSet` and whose value type is `GeomTraits::Vector_3`}
      \cgalParamDefault{`NormalMap()`}
   \cgalParamNBegin{k_neighbors}
     \cgalParamDescription{Shape detection: the number of neighbors for each point considered during region growing}
      \cgalParamType{`std::size_t`}
      \cgalParamDefault{12}
    \cgalParamNEnd
    \cgalParamNBegin{maximum_distance}
      \cgalParamDescription{Shape detection: the maximum distance from a point to a plane}
      \cgalParamType{`GeomTraits::FT`}
      \cgalParamDefault{2% of bounding box diagonal}
    \cgalParamNEnd
    \cgalParamNBegin{maximum_angle}
      \cgalParamDescription{Shape detection: maximum angle in degrees between the normal of a point and the plane normal}
      \cgalParamType{`GeomTraits::FT`}
      \cgalParamDefault{15 degrees}
    \cgalParamNEnd
    \cgalParamNBegin{minimum_region_size}
      \cgalParamDescription{Shape detection: minimum number of 3D points a region must have}
      \cgalParamType{`std::size_t`}
      \cgalParamDefault{1% of input points}
    \cgalParamNEnd
    \cgalParamNBegin{angle_tolerance}
      \cgalParamDescription{Shape regularization: maximum allowed angle in degrees between plane normals used for parallelism, orthogonality, and axis symmetry}
      \cgalParamType{`GeomTraits::FT`}
      \cgalParamDefault{5 degrees}
    \cgalParamNEnd
    \cgalParamNBegin{maximum_offset}
      \cgalParamDescription{Shape regularization: maximum allowed orthogonal distance between two parallel planes such that they are considered to be coplanar}
      \cgalParamType{`GeomTraits::FT`}
      \cgalParamDefault{0.5% of bounding box diagonal}
    \cgalParamNEnd
    \cgalParamNBegin{regularize_parallelism}
      \cgalParamDescription{Shape regularization: indicates whether parallelism should be regularized or not}
      \cgalParamType{boolean}
      \cgalParamDefault{false}
    \cgalParamNEnd
    \cgalParamNBegin{regularize_orthogonality}
      \cgalParamDescription{Shape regularization: indicates whether orthogonality should be regularized or not}
      \cgalParamType{boolean}
      \cgalParamDefault{false}
    \cgalParamNEnd
    \cgalParamNBegin{regularize_coplanarity}
      \cgalParamDescription{Shape regularization: indicates whether coplanarity should be regularized or not}
      \cgalParamType{boolean}
      \cgalParamDefault{true}
    \cgalParamNEnd
    \cgalParamNBegin{regularize_axis_symmetry}
      \cgalParamDescription{Shape regularization: indicates whether axis symmetry should be regularized or not}
      \cgalParamType{boolean}
      \cgalParamDefault{false}
    \cgalParamNEnd
    \cgalParamNBegin{symmetry_direction}
      \cgalParamDescription{Shape regularization: an axis for symmetry regularization}
      \cgalParamType{`GeomTraits::Vector_3`}
      \cgalParamDefault{Z axis that is `GeomTraits::Vector_3(0, 0, 1)`}
    \cgalParamNEnd
  \cgalNamedParamsEnd

  */
  template<typename CGAL_NP_TEMPLATE_PARAMETERS>
  std::size_t detect_planar_shapes(const CGAL_NP_CLASS& np = parameters::default_values()) {
    m_verbose = parameters::choose_parameter(parameters::get_parameter(np, internal_np::verbose), m_verbose);
    m_debug = parameters::choose_parameter(parameters::get_parameter(np, internal_np::debug), m_debug);

    if (m_verbose)
      std::cout << std::endl << "--- DETECTING PLANAR SHAPES: " << std::endl;

    m_planes.clear();
    m_polygons.clear();
    m_region_map.clear();

    m_point_map = Point_set_processing_3_np_helper<Point_range, CGAL_NP_CLASS, Point_map>::get_point_map(m_points, np);
    m_normal_map = Point_set_processing_3_np_helper<Point_range, CGAL_NP_CLASS, Normal_map>::get_normal_map(m_points, np);

    create_planar_shapes(np);

    CGAL_assertion(m_planes.size() == m_polygons.size());
    CGAL_assertion(m_polygons.size() == m_region_map.size());

    return m_polygons.size();
  }

  /*!
  \brief Retrieves the support planes of the detected and regularized shapes.

  @return
  vector with a `Plane_3` for each detected planar shape.

  \pre `shape detection performed`
  */
  const std::vector<Plane_3>& detected_planar_shapes() {
    return m_planes;
  }

  /*!
  \brief Retrieves the indices of detected and regularized shapes.

  @return
  indices into `points` for each detected planar shape.

  \pre `shape detection performed`
  */
  const std::vector<std::vector<std::size_t> >& detected_planar_shape_indices() {
    return m_planar_regions;
  }

  /*!
    \brief Detects and regularizes shapes in the provided point cloud and creates the kinetic space partition.

    Combines calls of `detect_planar_shapes()`, `initialize_partition()` and `partition()`.

    \tparam NamedParameters
    a sequence of \ref bgl_namedparameters "Named Parameters"

    \param k
    maximum number of allowed intersections for each input polygon before its expansion stops.

    \param np
    an instance of `NamedParameters`.

  \cgalNamedParamsBegin
    \cgalParamNBegin{point_map}
      \cgalParamDescription{a property map associating points to the elements of the point set `points`}
      \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type of the iterator of `PointSet` and whose value type is `GeomTraits::Point_3`}
          \cgalParamDefault{`PointMap()`}
    \cgalParamNEnd
    \cgalParamNBegin{normal_map}
      \cgalParamDescription{a property map associating normals to the elements of the point set `points`}
      \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type of the iterator of `PointSet` and whose value type is `GeomTraits::Vector_3`}
      \cgalParamDefault{`NormalMap()`}
   \cgalParamNBegin{k_neighbors}
     \cgalParamDescription{Shape detection: the number of neighbors for each point considered during region growing}
      \cgalParamType{`std::size_t`}
      \cgalParamDefault{12}
    \cgalParamNEnd
    \cgalParamNBegin{maximum_distance}
      \cgalParamDescription{Shape detection: the maximum distance from a point to a plane}
      \cgalParamType{`GeomTraits::FT`}
      \cgalParamDefault{2% of bounding box diagonal}
    \cgalParamNEnd
    \cgalParamNBegin{maximum_angle}
      \cgalParamDescription{Shape detection: maximum angle in degrees between the normal of a point and the plane normal}
      \cgalParamType{`GeomTraits::FT`}
      \cgalParamDefault{15 degrees}
    \cgalParamNEnd
    \cgalParamNBegin{minimum_region_size}
      \cgalParamDescription{Shape detection: minimum number of 3D points a region must have}
      \cgalParamType{`std::size_t`}
      \cgalParamDefault{1% of input points}
    \cgalParamNEnd
    \cgalParamNBegin{angle_tolerance}
      \cgalParamDescription{Shape regularization: maximum allowed angle in degrees between plane normals used for parallelism, orthogonality, and axis symmetry}
      \cgalParamType{`GeomTraits::FT`}
      \cgalParamDefault{5 degrees}
    \cgalParamNEnd
    \cgalParamNBegin{maximum_offset}
      \cgalParamDescription{Shape regularization: maximum allowed orthogonal distance between two parallel planes such that they are considered to be coplanar}
      \cgalParamType{`GeomTraits::FT`}
      \cgalParamDefault{0.5% of bounding box diagonal}
    \cgalParamNEnd
    \cgalParamNBegin{regularize_parallelism}
      \cgalParamDescription{Shape regularization: indicates whether parallelism should be regularized or not}
      \cgalParamType{boolean}
      \cgalParamDefault{false}
    \cgalParamNEnd
    \cgalParamNBegin{regularize_orthogonality}
      \cgalParamDescription{Shape regularization: indicates whether orthogonality should be regularized or not}
      \cgalParamType{boolean}
      \cgalParamDefault{false}
    \cgalParamNEnd
    \cgalParamNBegin{regularize_coplanarity}
      \cgalParamDescription{Shape regularization: indicates whether coplanarity should be regularized or not}
      \cgalParamType{boolean}
      \cgalParamDefault{true}
    \cgalParamNEnd
    \cgalParamNBegin{regularize_axis_symmetry}
      \cgalParamDescription{Shape regularization: indicates whether axis symmetry should be regularized or not}
      \cgalParamType{boolean}
      \cgalParamDefault{false}
    \cgalParamNEnd
    \cgalParamNBegin{symmetry_direction}
      \cgalParamDescription{Shape regularization: an axis for symmetry regularization}
      \cgalParamType{`GeomTraits::Vector_3`}
      \cgalParamDefault{Z axis that is `GeomTraits::Vector_3(0, 0, 1)`}
    \cgalParamNEnd
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
  \cgalNamedParamsEnd

  */
  template<typename CGAL_NP_TEMPLATE_PARAMETERS>
  void detection_and_partition(std::size_t k, const CGAL_NP_CLASS& np = parameters::default_values()) {
    detect_planar_shapes(np);
    initialize_partition(np);
    partition(k);
  }

  /*!
  \brief initializes the kinetic partition.

  \param np
  a sequence of \ref bgl_namedparameters "Named Parameters"
  among the ones listed below

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
  \cgalNamedParamsEnd

    \pre shape detection performed
  */
  template<typename CGAL_NP_TEMPLATE_PARAMETERS>
  void initialize_partition(const CGAL_NP_CLASS& np = parameters::default_values()) {
    m_kinetic_partition.insert(m_polygon_pts, m_polygon_indices, np);

    m_kinetic_partition.initialize(np);
  }

  /*!
  \brief Propagates the kinetic polygons in the initialized partition.

  \param k
  maximum number of allowed intersections for each input polygon before its expansion stops.

  \pre partition initialized
  */
  void partition(std::size_t k) {
    FT partition_time, finalization_time, conformal_time;
    m_kinetic_partition.partition(k, partition_time, finalization_time, conformal_time);
    std::cout << "Bounding box partitioned into " << m_kinetic_partition.number_of_volumes() << " volumes" << std::endl;

    m_kinetic_partition.get_linear_cell_complex(m_lcc);

    setup_energyterms();
  }

  /*!
  \brief Access to the kinetic partition.

  @return
  created kinetic partition data structure

  \pre partition created
  */
  const KSP& kinetic_partition() const {
    return m_kinetic_partition;
  }

  /*!
  \brief Uses min-cut to solve an inside/outside labeling of the volumes of the kinetic partition and provides the reconstructed surface as a list of indexed polygons.
  Estimates a horizontal ground plane within the detected shapes. Cells in the partition below the ground plane receive a weight to be labeled as inside.
  The z axis is considered as vertical upwards pointing.

  \tparam OutputPointIterator
  an output iterator taking Point_3.

  \tparam OutputPolygonIterator
  an output iterator taking polygon indices std::vector<std::size_t>.

  \param lambda
  trades data faithfulness of the reconstruction for low complexity. Must be in the range `[0, 1)`.

  \param pit
  output iterator to receive the vertices of the reconstructed surface.

  \param polyit
  output iterator to store all polygonal faces of the reconstructed surface.

  \pre partition created
  */
  template<class OutputPointIterator, class OutputPolygonIterator>
  void reconstruct_with_ground(FT lambda, OutputPointIterator pit, OutputPolygonIterator polyit) {
    KSR_3::Graphcut<Kernel> gc(lambda);

    // add ground consideration here
    // m_cost_matrix and m_face_neighbors_lcc should contain the original values without consideration of ground/preset cell labels
    // Create copy here if necessary and fill both vectors according to parameters of reconstruction

    set_outside_volumes(true, m_cost_matrix);

    std::vector<std::pair<std::size_t, std::size_t> > edges(m_face_neighbors_lcc.size());

    std::size_t redirected = 0;
    for (std::size_t i = 0; i < m_face_neighbors_lcc.size(); i++) {
      std::size_t lower, upper;
      if (m_face_neighbors_lcc[i].first < m_face_neighbors_lcc[i].second) {
        lower = m_face_neighbors_lcc[i].first;
        upper = m_face_neighbors_lcc[i].second - 6;
      }
      else {
        lower = m_face_neighbors_lcc[i].second;
        upper = m_face_neighbors_lcc[i].first - 6;
      }
      // Check if the face is on a bbox face besides the top face.
      // If so and the connected volume is below the ground, redirect the face to the bottom face node.
      if ((lower < 6) && m_volume_below_ground[upper]) {
        redirected++;

        if (m_face_neighbors_lcc[i].second < 6) {
          edges[i].first = m_face_neighbors_lcc[i].first;
          edges[i].second = 0;
        }

        if (m_face_neighbors_lcc[i].first < 6) {
          m_face_neighbors_lcc[i].first = 0;
          edges[i].second = m_face_neighbors_lcc[i].second;
        }
      }
      else edges[i] = m_face_neighbors_lcc[i];
    }

    if (m_verbose)
      std::cout << redirected << " faces redirected to below ground" << std::endl;

    gc.solve(edges, m_face_area_lcc, m_cost_matrix, m_labels);

    reconstructed_model_polylist_lcc(pit, polyit, lambda);
  }

  /*!
  \brief Uses min-cut to solve an inside/outside labeling of the volumes of the kinetic partition and provides the reconstructed surface as a list of indexed polygons.
  The `external_nodes` parameter allows to indicate the preferred labels for faces on the bounding box.

  \tparam OutputPointIterator
  an output iterator taking `Point_3`.

  \tparam OutputPolygonIterator
  an output iterator taking polygon indices `std::vector<std::size_t>`.

  \param lambda
  trades data faithfulness of the reconstruction for low complexity. Should be in the range [0, 1).

  \param external_nodes
  adds label preference for the faces on the bounding box. Bounding box sides without preset label are chosen by the min-cut.
  Setting `external_nodes[ZMIN] = true` sets the inside label as the preferred label for the ZMIN side of the bounding box.

  \param pit
  output iterator to receive the vertices of the reconstructed surface.

  \param polyit
  output iterator to store all polygonal faces of the reconstructed surface.

  \pre partition created
  */
  template<class OutputPointIterator, class OutputPolygonIterator>
  void reconstruct(FT lambda, std::map<typename KSP::Face_support, bool> external_nodes, OutputPointIterator pit, OutputPolygonIterator polyit) {
    KSR_3::Graphcut<Kernel> gc(lambda);

    // add node consideration here
    set_outside_volumes(false, m_cost_matrix);

    const std::size_t force = m_total_inliers * 3;
    for (auto& p : external_nodes) {
      int idx = -p.first - 1;
      if (p.second) {
        m_cost_matrix[0][idx] = force;
        m_cost_matrix[1][idx] = 0;
      }
      else {
        m_cost_matrix[0][idx] = 0;
        m_cost_matrix[1][idx] = force;
      }
    }

    gc.solve(m_face_neighbors_lcc, m_face_area_lcc, m_cost_matrix, m_labels);

    reconstructed_model_polylist_lcc(pit, polyit, lambda);
  }

private:
  using Point_2 = typename Kernel::Point_2;
  using Vector_3 = typename Kernel::Vector_3;
  using Triangle_2 = typename Kernel::Triangle_2;

  using Indices = std::vector<std::size_t>;
  using Polygon_3 = std::vector<Point_3>;

  using Region_type = CGAL::Shape_detection::Point_set::Least_squares_plane_fit_region_for_point_set<Point_range>;
  using Neighbor_query = CGAL::Shape_detection::Point_set::K_neighbor_query_for_point_set<Point_range>;
  using Sorting = CGAL::Shape_detection::Point_set::Least_squares_plane_fit_sorting_for_point_set<Point_range, Neighbor_query>;
  using Region_growing = CGAL::Shape_detection::Region_growing<Neighbor_query, Region_type>;
  using From_exact = typename CGAL::Cartesian_converter<Intersection_kernel, Kernel>;

  struct Vertex_info { FT z = FT(0); };
  struct Face_info { };

  using Fbi = CGAL::Triangulation_face_base_with_info_2<Face_info, Kernel>;
  //using Fb = CGAL::Alpha_shape_face_base_2<Kernel, Fbi>;

  using Vbi = CGAL::Triangulation_vertex_base_with_info_2<Vertex_info, Kernel>;
  //using Vb = CGAL::Alpha_shape_vertex_base_2<Kernel, Vbi>;

  using Tds = CGAL::Triangulation_data_structure_2<Vbi, Fbi>;
  using Delaunay_2 = CGAL::Delaunay_triangulation_2<Kernel, Tds>;

  using Delaunay_3 = CGAL::Delaunay_triangulation_3<Kernel>;

  typedef CGAL::Linear_cell_complex_traits<3, CGAL::Exact_predicates_exact_constructions_kernel> Traits;
  using LCC = CGAL::Linear_cell_complex_for_combinatorial_map<3, 3, Traits, typename KSP::Linear_cell_complex_min_items>;
  using Dart_descriptor = typename LCC::Dart_descriptor;
  using Dart = typename LCC::Dart;

  struct VI {
    VI() : i(-1), j(-1) {}
    int i, j;
    Dart_descriptor dh;
    typename Intersection_kernel::Point_2 p;
  };

  typedef CGAL::Triangulation_vertex_base_with_info_2<VI, Intersection_kernel> Vbi2;
  typedef CGAL::Constrained_triangulation_face_base_2<Intersection_kernel>   Fb;
  typedef CGAL::Triangulation_data_structure_2<Vbi2, Fb>       Tds2;
  typedef CGAL::Exact_intersections_tag  Itag;
  typedef CGAL::Constrained_Delaunay_triangulation_2<Intersection_kernel, Tds2, Itag> CDT;

  typedef typename CDT::Vertex_handle            Vertex_handle;
  typedef typename CDT::Face_handle              Face_handle;
  typedef typename CDT::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename CDT::Finite_edges_iterator    Finite_edges_iterator;
  typedef typename CDT::Finite_faces_iterator    Finite_faces_iterator;

  //using Visibility = KSR_3::Visibility<Kernel, Intersection_kernel, Point_map, Normal_map>;
  using Index = typename KSP::Index;
  using Face_attribute = typename LCC::Base::template Attribute_descriptor<2>::type;
  using Volume_attribute = typename LCC::Base::template Attribute_descriptor<3>::type;

  bool m_verbose;
  bool m_debug;

  Point_range &m_points;
  Point_map m_point_map;
  Normal_map m_normal_map;

  std::vector<std::vector<std::size_t> > m_planar_regions;
  std::vector<typename Region_growing::Primitive_and_region> m_regions;
  std::map<std::size_t, Indices> m_region_map;
  double m_detection_distance_tolerance;

  std::size_t m_ground_polygon_index;
  Plane_3 m_ground_plane;

  std::vector<Plane_3> m_planes;
  std::vector<Point_3> m_polygon_pts;
  std::vector<std::vector<std::size_t> > m_polygon_indices;
  std::vector<Polygon_3> m_polygons;
  KSP m_kinetic_partition;

  LCC m_lcc;
  std::vector<typename LCC::Dart_const_descriptor> m_faces_lcc;
  std::map<Face_attribute, std::size_t> m_attrib2index_lcc;
  std::vector<std::size_t> lcc2index;
  std::vector<std::size_t> index2lcc;

  // Face indices are now of type Indices and are not in a range 0 to n
  std::vector<Indices> m_face_inliers;
  std::vector<FT> m_face_area, m_face_area_lcc;
  std::vector<std::pair<std::size_t, std::size_t> > m_face_neighbors_lcc;
  std::map<std::pair<std::size_t, std::size_t>, std::size_t> m_neighbors2index_lcc;

  std::vector<std::pair<std::size_t, std::size_t> > m_volume_votes; // pair<inside, outside> votes
  std::vector<bool> m_volume_below_ground;
  std::vector<std::vector<double> > m_cost_matrix;
  std::vector<FT> m_volumes; // normalized volume of each kinetic volume
  std::vector<std::size_t> m_labels;

  std::size_t m_total_inliers;


  /*!
  \brief Creates the visibility (data-) and regularity energy terms from the input point cloud and the kinetic partition.

  \pre `successful initialization`
  */
  void setup_energyterms() {
    if (m_lcc.template one_dart_per_cell<3>().size() == 0) {
      std::cout << "Kinetic partition is not constructed or does not have volumes" << std::endl;
      return;
    }

    m_face_area.clear();
    m_face_inliers.clear();

    auto face_range = m_lcc.template one_dart_per_cell<2>();
    m_faces_lcc.reserve(face_range.size());

    for (auto& d : face_range) {
      typename LCC::Dart_descriptor dh = m_lcc.dart_descriptor(d);

      Face_attribute fa = m_lcc.template attribute<2>(dh);
      if (fa == m_lcc.null_descriptor) {
        dh = m_lcc.template beta<3>(dh);
        fa = m_lcc.template attribute<2>(dh);
      }

      if (fa == m_lcc.null_descriptor) {
        std::cout << "null dart 1 " << m_lcc.template one_dart_per_incident_cell<3, 2>(dh).size() << std::endl;
      }

      m_faces_lcc.push_back(dh);

      auto p = m_attrib2index_lcc.emplace(std::make_pair(m_lcc.template attribute<2>(m_faces_lcc.back()), m_faces_lcc.size() - 1));
      CGAL_assertion(p.second);
    }

    // Create value arrays for graph cut
    m_face_inliers.resize(m_faces_lcc.size());
    m_face_area.resize(m_faces_lcc.size());
    m_face_area_lcc.resize(m_faces_lcc.size());
    m_face_neighbors_lcc.resize(m_faces_lcc.size(), std::pair<int, int>(-1, -1));

    m_cost_matrix.resize(2);
    m_cost_matrix[0].resize(m_kinetic_partition.number_of_volumes() + 6, 0);
    m_cost_matrix[1].resize(m_kinetic_partition.number_of_volumes() + 6, 0);

    for (std::size_t i = 0; i < m_faces_lcc.size(); i++) {
      auto n = m_lcc.template one_dart_per_incident_cell<3, 2>(m_faces_lcc[i]);

      assert(n.size() == 1 || n.size() == 2);
      auto it = n.begin();

      auto& finf = m_lcc.template info<2>(m_faces_lcc[i]);

      bool skipped = false;

      Volume_attribute va = m_lcc.template attribute<3>(m_lcc.dart_descriptor(*it));
      if (va == m_lcc.null_descriptor) {
        skipped = true;
        it++;
      }

      if (it == n.end()) {
        std::cout << "face not connected to a volume" << std::endl;
        continue;
      }

      va = m_lcc.template attribute<3>(m_lcc.dart_descriptor(*it));
      if (va == m_lcc.null_descriptor) {
        write_face(m_lcc.dart_descriptor(*it), "face_wo_volume.ply");
      }

      int first = static_cast<int>(m_lcc.template info<3>(m_lcc.dart_descriptor(*it)).volume_id);
      auto& inf1 = m_lcc.template info<3>(m_lcc.dart_descriptor(*it++));

      auto inf2 = inf1;
      if (n.size() == 2 && it != n.end())
        inf2 = m_lcc.template info<3>(m_lcc.dart_descriptor(*it));

      int second;
      if (n.size() == 2 && it != n.end())
        second = static_cast<int>(m_lcc.template info<3>(m_lcc.dart_descriptor(*it)).volume_id);

      if (n.size() == 2 && it != n.end())
        m_face_neighbors_lcc[i] = std::make_pair(first + 6, m_lcc.template info<3>(m_lcc.dart_descriptor(*it)).volume_id + 6);
      else
        m_face_neighbors_lcc[i] = std::make_pair(first + 6, -m_lcc.template info<2>(m_faces_lcc[i]).input_polygon_index - 1);

      if (m_face_neighbors_lcc[i].first > m_face_neighbors_lcc[i].second)
        m_face_neighbors_lcc[i] = std::make_pair(m_face_neighbors_lcc[i].second, m_face_neighbors_lcc[i].first);

      if (m_face_neighbors_lcc[i].first < m_face_neighbors_lcc[i].second) {
        auto it = m_neighbors2index_lcc.emplace(std::make_pair(m_face_neighbors_lcc[i], i));
        assert(it.second);
      }
    }

    check_ground();

    m_face_inliers.clear();
    m_face_inliers.resize(m_faces_lcc.size());
    collect_points_for_faces_lcc();
    count_volume_votes_lcc();

    if (m_verbose)
      std::cout << "* computing data term ... ";

    std::size_t max_inside = 0, max_outside = 0;
    for (std::size_t i = 0; i < m_volumes.size(); i++) {
      max_inside = (std::max<std::size_t>)(static_cast<std::size_t>(m_cost_matrix[0][i + 6]), max_inside);
      max_outside = (std::max<std::size_t>)(static_cast<std::size_t>(m_cost_matrix[1][i + 6]), max_outside);
    }

    // Dump volumes colored by votes
/*
    if (false) {
      namespace fs = boost::filesystem;
      for (fs::directory_iterator end_dir_it, it("gc/i"); it != end_dir_it; ++it) {
        fs::remove_all(it->path());
      }
      for (fs::directory_iterator end_dir_it, it("gc/o"); it != end_dir_it; ++it) {
        fs::remove_all(it->path());
      }
      for (fs::directory_iterator end_dir_it, it("gc/n"); it != end_dir_it; ++it) {
        fs::remove_all(it->path());
      }
      for (fs::directory_iterator end_dir_it, it("gc/all"); it != end_dir_it; ++it) {
        fs::remove_all(it->path());
      }
      for (std::size_t i = 0; i < m_volumes.size(); i++) {
        // skip 0/0 volumes? Maybe safe them a few seconds later to be able to separate them?
        CGAL::Color c;

        int diff = int(m_cost_matrix[0][i + 6]) - int(m_cost_matrix[1][i + 6]);

        if (diff > 0) {
          std::size_t m = (std::max<int>)(50, (diff * 255) / max_inside);
          c = CGAL::Color(0, m, 0);
        }
        else {
          std::size_t m = (std::max<int>)(50, (-diff * 255) / max_outside);
          c = CGAL::Color(0, 0, m);
        }

        if (diff < 0) {
          dump_volume(i, "gc/o/" + std::to_string(i) + "-vol-" + std::to_string(m_cost_matrix[0][i + 6]) + "-" + std::to_string(m_cost_matrix[1][i + 6]), c);
          dump_volume(i, "gc/all/" + std::to_string(i) + "-vol-" + std::to_string(m_cost_matrix[0][i + 6]) + "-" + std::to_string(m_cost_matrix[1][i + 6]), c);
        }
        else if (diff > 0) {
          dump_volume(i, "gc/i/" + std::to_string(i) + "-vol-" + std::to_string(m_cost_matrix[0][i + 6]) + "-" + std::to_string(m_cost_matrix[1][i + 6]), c);
          dump_volume(i, "gc/all/" + std::to_string(i) + "-vol-" + std::to_string(m_cost_matrix[0][i + 6]) + "-" + std::to_string(m_cost_matrix[1][i + 6]), c);
        }
        else {
          dump_volume(i, "gc/n/" + std::to_string(i) + "-vol-0-0", CGAL::Color(255, 255, 255));
          dump_volume(i, "gc/all/" + std::to_string(i) + "-vol-0-0", CGAL::Color(255, 255, 255));
        }
      }
    }
*/
  }

  /*!
  \brief Provides the data and regularity energy terms for reconstruction via min-cut.

  \param edges
  contains a vector of pairs of volume indices. Indicates which volumes should be connected in the min-cut formulation.

  \param edge_costs
  contains the cost for each edge specified in `edges` for two labels with different labels. For equal labels, the cost is 0. Needs to be index compatible to the `edges` parameter.

  \param cost_matrix
  provides the cost of a label for each volume cell. The first index corresponds to the label and the second index corresponds to the volume index.

  @return
  fails if the dimensions of parameters does not match the kinetic partition.

  \pre `successful initialization`
  */
  template<typename NamedParameters>
  bool setup_energyterms(
    const std::vector< std::pair<std::size_t, std::size_t> >& edges,
    const std::vector<double>& edge_costs,
    const std::vector< std::vector<double> >& cost_matrix);

  /*!
  \brief Provides the reconstructed surface as a list of indexed triangles.

  \param pit
  an output iterator taking Point_3.

  \param triit
  an output iterator taking std::size_t.

  \pre `successful reconstruction`
  */
  template<class OutputPointIterator, class OutputTriangleIterator>
  void reconstructed_model_trilist(OutputPointIterator pit, OutputTriangleIterator triit) {
    if (m_labels.empty())
      return;

    std::map<Point_3, std::size_t> pt2idx;

    for (std::size_t i = 0; i < m_faces_lcc.size(); i++) {
      const auto& n = m_face_neighbors_lcc[i];
      // Do not extract boundary faces.
      if (n.second < 6)
        continue;
      if (m_labels[n.first] != m_labels[n.second]) {
        std::vector<Point_3> face;
        m_kinetic_partition.vertices(m_faces_lcc[i], std::back_inserter(face));

        std::vector<std::size_t> indices(face.size());

        for (std::size_t i = 0; i < face.size(); i++) {
          auto p = pt2idx.emplace(face[i], pt2idx.size());
          if (p.second)
            *pit++ = face[i];
          indices[i] = p.first->second;
        }

        for (std::size_t i = 2; i < face.size(); i++) {
          *triit++ = indices[0];
          *triit++ = indices[i - 1];
          *triit++ = indices[i];
        }
      }
    }
  }

  /*!
  \brief Provides the reconstructed surface as a list of indexed polygons.

  \param pit
  an output iterator taking Point_3.

  \param triit
  an output iterator taking std::vector<std::size_t>.

  \pre `successful reconstruction`
  */
  template<class OutputPointIterator, class OutputPolygonIterator>
  void reconstructed_model_polylist_lcc(OutputPointIterator pit, OutputPolygonIterator polyit, FT lambda) {
    if (m_labels.empty())
      return;

    From_exact from_exact;

    std::map<Point_3, std::size_t> pt2idx;

    std::vector<int> region_index(m_faces_lcc.size(), -1);
    std::size_t region = 0;

    std::vector<std::vector<std::vector<Point_3> > > polygon_regions;
    std::vector<typename Intersection_kernel::Plane_3> planes;

    for (std::size_t i = 0; i < m_faces_lcc.size(); i++) {
      const auto& n = m_face_neighbors_lcc[i];
      //       if (n.first < 6 || n.second < 6)
      //         continue;

      if (m_labels[n.first] != m_labels[n.second]) {
        Face_attribute fa = m_lcc.template attribute<2>(m_faces_lcc[i]);

        if (fa == m_lcc.null_descriptor)
          std::cout << "null dart 1" << std::endl;

        if (m_labels[m_lcc.template info<3>(m_faces_lcc[i]).volume_id + 6] == 0) {
          Dart_descriptor dh = m_lcc.template beta<3>(m_faces_lcc[i]);
          if (dh == m_lcc.null_dart_descriptor)
            continue;
          if (m_lcc.template attribute<3>(m_faces_lcc[i]) == m_lcc.null_descriptor)
            continue;
          m_faces_lcc[i] = dh;
        }

        fa = m_lcc.template attribute<2>(m_faces_lcc[i]);

        if (fa == m_lcc.null_descriptor) {
          std::cout << "null dart 2" << std::endl;
          continue;
        }

        if (region_index[fa] == -1) {
          std::vector<std::vector<Point_3> > faces;

          collect_connected_component(m_faces_lcc[i], region_index, region++, faces);
          planes.push_back(m_lcc.template info_of_attribute<2>(fa).plane);
          polygon_regions.push_back(std::move(faces));
        }
      }
    }

    for (std::size_t i = 0; i < m_faces_lcc.size(); i++) {
      Face_attribute fa = m_lcc.template attribute<2>(m_faces_lcc[i]);
      if (region_index[fa] == -1)
        continue;

      if (m_labels[m_lcc.template info<3>(m_faces_lcc[i]).volume_id + 6] == 0)
        std::cout << "outside face" << std::endl;
    }

    if (m_verbose)
      std::cout << "polygon regions " << polygon_regions.size() << std::endl;

    if (m_debug)
      KSP_3::internal::dump_polygons(polygon_regions, "faces_by_region-" + std::to_string(lambda) + ".ply");

    std::vector<std::vector<std::size_t> > borders;
    std::vector<std::vector<std::size_t> > borders_per_region;
    collect_connected_border(borders, region_index, borders_per_region);

    for (std::size_t i = 0; i < region; i++) {
      if (borders_per_region[i].size() > 0) {
        /*std::size_t outer = -1;
        typename Intersection_kernel::FT min = (std::numeric_limits<double>::max)();
        for (std::size_t j = 0; j < borders_per_region[i].size(); j++)
          for (std::size_t k = 0; k < borders[borders_per_region[i][j]].size(); k++) {
            const typename Intersection_kernel::Point_3& p = m_lcc.point(m_lcc.dart_of_attribute<0>(borders[borders_per_region[i][j]][k]));
            if (p.x() < min) {
              min = p.x();
              outer = j;
            }
          }

        for (std::size_t j = 0; j < borders_per_region[i].size(); j++) {
          std::string fn;
          if (j == outer)
            fn = std::to_string(i) + "-outer.polylines.txt";
          else
            fn = std::to_string(i) + "-" + std::to_string(j) + ".polylines.txt";
          std::ofstream vout(fn);
          vout << (borders[borders_per_region[i][j]].size() + 1);
          for (std::size_t k = 0; k < borders[borders_per_region[i][j]].size(); k++) {
            vout << " " << from_exact(m_lcc.point(m_lcc.dart_of_attribute<0>(borders[borders_per_region[i][j]][k])));
          }
          vout << " " << from_exact(m_lcc.point(m_lcc.dart_of_attribute<0>(borders[borders_per_region[i][j]][0]))) << std::endl;
          vout.close();
        }*/

        if (borders_per_region[i].size() > 1) {
          std::vector<std::vector<std::size_t> > polygons;
          polygons.reserve(borders_per_region[i].size());
          for (std::size_t j = 0; j < borders_per_region[i].size(); j++)
            polygons.push_back(std::move(borders[borders_per_region[i][j]]));

          insert_ghost_edges_cdt(polygons, planes[i]);

          borders_per_region[i].resize(1);
          CGAL_assertion(borders[borders_per_region[i][0]].empty());
          CGAL_assertion(!polygons[0].empty());
          borders[borders_per_region[i][0]] = std::move(polygons[0]);
        }
      }
    }

    std::map<std::size_t, std::size_t> attrib2idx;
    for (std::size_t i = 0; i < borders.size(); i++) {
      if (borders[i].empty())
        continue;

      std::vector<std::size_t> indices(borders[i].size());
      for (std::size_t j = 0; j != borders[i].size(); j++) {
        auto p = attrib2idx.emplace(borders[i][j], attrib2idx.size());
        if (p.second)
          *pit++ = from_exact(m_lcc.point(m_lcc.template dart_of_attribute<0>(borders[i][j])));
        indices[j] = p.first->second;
      }

      std::reverse(indices.begin(), indices.end());

      *polyit++ = std::move(indices);
    }
  }

  /*!
  \brief Provides the reconstructed surface as a list of indexed polygons.

  \param pit
  an output iterator taking Point_3.

  \param triit
  an output iterator taking std::vector<std::size_t>.

  \pre `successful reconstruction`
  */
  template<class OutputPointIterator, class OutputPolygonIterator>
  void reconstructed_model_polylist(OutputPointIterator pit, OutputPolygonIterator polyit) {
    if (m_labels.empty())
      return;

    From_exact from_exact;

    std::map<Point_3, std::size_t> pt2idx;

    for (std::size_t i = 0; i < m_faces_lcc.size(); i++) {
      const auto& n = m_face_neighbors_lcc[i];
      // Do not extract boundary faces.
      if (n.second < 6)
        continue;
      if (m_labels[n.first] != m_labels[n.second]) {
        std::vector<Point_3> face;

        for (const auto& vd : m_lcc.one_dart_per_incident_cell<0, 2>(m_faces_lcc[i]))
          face.push_back(from_exact(m_lcc.point(m_lcc.dart_descriptor(vd))));

        std::vector<std::size_t> indices(face.size());

        for (std::size_t i = 0; i < face.size(); i++) {
          auto p = pt2idx.emplace(face[i], pt2idx.size());
          if (p.second)
            *pit++ = face[i];
          indices[i] = p.first->second;
        }

        *polyit++ = std::move(indices);
      }
    }
  }

  std::size_t add_convex_hull_shape(
    const std::vector<std::size_t>& region, const Plane_3& plane) {

    std::vector<Point_2> points;
    points.reserve(region.size());
    for (const std::size_t idx : region) {
      CGAL_assertion(idx < m_points.size());
      const auto& p = get(m_point_map, idx);
      const auto q = plane.projection(p);
      const auto point = plane.to_2d(q);
      points.push_back(point);
    }
    CGAL_assertion(points.size() == region.size());

    std::vector<Point_2> ch;
    CGAL::convex_hull_2(points.begin(), points.end(), std::back_inserter(ch));

    std::vector<Point_3> polygon;
    for (const auto& p : ch) {
      const auto point = plane.to_3d(p);
      polygon.push_back(point);
    }

    const std::size_t shape_idx = m_polygons.size();
    m_polygons.push_back(polygon);
    m_planes.push_back(plane);

    m_polygon_indices.push_back(std::vector<std::size_t>());
    m_polygon_indices.back().resize(polygon.size());
    std::iota(std::begin(m_polygon_indices.back()), std::end(m_polygon_indices.back()), m_polygon_pts.size());
    std::copy(polygon.begin(), polygon.end(), std::back_inserter(m_polygon_pts));

    return shape_idx;
  }

  void store_convex_hull_shape(const std::vector<std::size_t>& region, const Plane_3& plane, std::vector<std::vector<Point_3> > &polys) {

    std::vector<Point_2> points;
    points.reserve(region.size());
    for (const std::size_t idx : region) {
      CGAL_assertion(idx < m_points.size());
      const auto& p = get(m_point_map, idx);
      const auto q = plane.projection(p);
      const auto point = plane.to_2d(q);
      points.push_back(point);
    }
    CGAL_assertion(points.size() == region.size());

    std::vector<Point_2> ch;
    CGAL::convex_hull_2(points.begin(), points.end(), std::back_inserter(ch));

    std::vector<Point_3> polygon;
    for (const auto& p : ch) {
      const auto point = plane.to_3d(p);
      polygon.push_back(point);
    }

    polys.push_back(polygon);
  }

  std::pair<int, int> make_canonical_pair(int i, int j)
  {
    if (i > j) return std::make_pair(j, i);
    return std::make_pair(i, j);
  }

  void check_ground() {
    std::size_t num_volumes = m_kinetic_partition.number_of_volumes();
    // Set all volumes to not be below the ground, this leads to the standard 6 outside node connection.
    m_volume_below_ground.resize(num_volumes, false);
    From_exact from_exact;

    if (m_ground_polygon_index != static_cast<std::size_t>(-1))
      for (const auto &vd : m_lcc.template one_dart_per_cell<3>()) {
        const auto& info = m_lcc. template info<3>(m_lcc.dart_descriptor(vd));

        m_volume_below_ground[info.volume_id] = (from_exact(info.barycenter) - m_regions[m_ground_polygon_index].first.projection(from_exact(info.barycenter))).z() < 0;
      }
  }

  void collect_connected_component(typename LCC::Dart_descriptor face, std::vector<int>& region_index, std::size_t region, std::vector<std::vector<Point_3> > &faces) {
    std::queue<std::size_t> face_queue;
    face_queue.push(face);

    From_exact from_exact;

    auto& finfo = m_lcc.template info<2>(face);
    int ip = m_lcc.template info<2>(face).input_polygon_index;
    typename Intersection_kernel::Plane_3 pl = m_lcc.template info<2>(face).plane;

    if (m_labels[m_lcc.template info<3>(face).volume_id + 6] == 0)
      std::cout << "collect_connected_component called on outside face" << std::endl;

    while (!face_queue.empty()) {
      Dart_descriptor cur_fdh(face_queue.front());
      Face_attribute cur_fa = m_lcc.template attribute<2>(cur_fdh);
      face_queue.pop();

      if (region_index[cur_fa] == region)
        continue;

      //write_face(cur_fdh, std::to_string(region) + "-inside-" + std::to_string(cur_fa) + ".ply");

      region_index[cur_fa] = static_cast<int>(region);

      Dart_descriptor n = cur_fdh;
      std::vector<Point_3> f;
      do {
        f.push_back(from_exact(m_lcc.point(n)));
        n = m_lcc.beta(n, 1);
      } while (n != cur_fdh);
      faces.push_back(std::move(f));

      // Iterate over edges of face.

      Dart_descriptor edh = cur_fdh;
      do {
        Dart_descriptor fdh = m_lcc.beta(edh, 2, 3);
        do {
          Face_attribute fa = m_lcc.template attribute<2>(fdh);

          if (fa == m_lcc.null_descriptor) {
            // fdh is outside of the bbox, switching back inside to check the face on the boundary
            fdh = m_lcc.template beta<3>(fdh);
            fa = m_lcc.template attribute<2>(fdh);

            if (fa == m_lcc.null_descriptor)
              break;
          }

          auto& finfo2 = m_lcc.template info<2>(fdh);
          if (fa == cur_fa) {
            fdh = m_lcc.template beta<2, 3>(fdh);
            continue;
          }
          auto& inf = m_lcc.template info<2>(fdh);
          bool added = false;

          //write_face(fdh, std::to_string(region) + "-" + std::to_string(fa) + ".ply");

          const auto& n = m_face_neighbors_lcc[m_attrib2index_lcc[fa]];

          // Belongs to reconstruction?
          bool internal = m_labels[n.first] == m_labels[n.second];
          if (m_labels[n.first] == m_labels[n.second]) {
            fdh = m_lcc.template beta<2>(fdh);
            Dart_descriptor fdh2 = m_lcc.template beta<3>(fdh);
            if (fdh2 != m_lcc.null_dart_descriptor)
              fdh = fdh2;

            continue;
          }

          // Already segmented?
          if (region_index[fa] != -1) {
            if (!internal)
              break;

            fdh = m_lcc.template beta<2>(fdh);
            Dart_descriptor fdh2 = m_lcc.template beta<3>(fdh);
            if (fdh2 != m_lcc.null_dart_descriptor)
              fdh = fdh2;

            continue;
          }

          // If the face is part of the reconstruction, but on the inside volume, switch to the mirror face on the outside.
          if (n.first >= 6 && n.second >= 6 && m_labels[m_lcc.template info<3>(fdh).volume_id + 6] == 0) {
            fdh = m_lcc.template beta<3>(fdh);
            fa = m_lcc.template attribute<2>(fdh);
            finfo2 = m_lcc.template info<2>(fdh);
          }

          if (ip != -7) {
            if (m_lcc.template info<2>(fdh).input_polygon_index == ip) {
              if (internal)
                break;

              added = true;
              face_queue.push(fdh);

//               if (debug)
//                 std::cout << ip << std::endl;
//
//               if (debug)
//                 write_face(fdh, std::to_string(region) + "-inside-" + std::to_string(fa) + ".ply");
            }
            else
              if (!internal)
                break;
          }
          else
            if (m_lcc.template info<2>(fdh).plane == pl || m_lcc.template info<2>(fdh).plane == pl.opposite()) {
              if (internal)
                break;

              added = true;
              Plane_3 pla = from_exact(pl);

//               if (debug)
//                 std::cout << ip << " " << pl.a() << " " << pl.b() << " " << pl.c() << " " << pl.d() << std::endl;
//
//               if (debug)
//                 write_face(fdh, std::to_string(region) + "-inside-" + std::to_string(fa) + ".ply");

              face_queue.push(fdh);
            }
            else
              if (!internal)
                break;

//           if (!added)
//             border_edges.push_back(edh);

          break;
        } while (fdh != edh);
        edh = m_lcc.template beta<1>(edh);
      } while (edh != cur_fdh);
    }
  }

  bool is_border_edge(typename LCC::Dart_descriptor dh) {
    const Face_attribute& fa = m_lcc.template attribute<2>(dh);
    auto& finfo = m_lcc.template info_of_attribute<2>(fa);

    if (!m_labels[m_lcc.template info<3>(dh).volume_id + 6] == 1) {
      write_face(dh, "flipface.ply");
      std::cout << "is_border_edge called on dart of outside volume, dh " << dh << " volume_id " << m_lcc.template info<3>(dh).volume_id << std::endl;
    }

    Dart_descriptor edh = m_lcc.beta(dh, 2, 3);
    do {
      Face_attribute fa2 = m_lcc.template attribute<2>(edh);
      if (fa2 == m_lcc.null_descriptor)
        return true;

//       if (debug)
//         write_face(edh, "cur_is_border.ply");

      if (fa2 == fa) {
        std::cout << "should not happen" << std::endl;
        edh = m_lcc.template beta<2, 3>(edh);
        continue;
      }

      const auto& n = m_face_neighbors_lcc[m_attrib2index_lcc[fa2]];
      bool internal = (m_labels[n.first] == m_labels[n.second]);

      auto& finfo2 = m_lcc.template info_of_attribute<2>(fa2);
      // Is neighbor face on same support plane?
      if (finfo2.input_polygon_index != finfo.input_polygon_index)
        if (!internal)
          return true;
        else {

          edh = m_lcc.template beta<2>(edh);
          Dart_descriptor edh2 = m_lcc.template beta<3>(edh);
          if (edh2 != m_lcc.null_dart_descriptor)
            edh = edh2;

          continue;
        }

      if (finfo2.input_polygon_index == -7)
        if (finfo2.plane != finfo.plane && finfo2.plane != finfo.plane.opposite())
          if (!internal)
            return true;
          else {
            edh = m_lcc.template beta<2>(edh);
            Dart_descriptor edh2 = m_lcc.template beta<3>(edh);
            if (edh2 != m_lcc.null_dart_descriptor)
              edh = edh2;

            continue;
          };
      return internal;
    } while (edh != dh);

    // If there is no neighbor face on the same support plane, this is a border edge.
    return true;
  }

  void insert_ghost_edges_cdt(std::vector<std::vector<std::size_t> >& polygons, const typename Intersection_kernel::Plane_3 pl) const {
    CDT cdt;
    From_exact from_exact;

    std::unordered_map<std::size_t, std::size_t> va2vh;
    std::vector<Vertex_handle> vertices;

    std::size_t num_vertices = 0;

    for (int i = 0; i < polygons.size(); i++) {
      num_vertices += polygons[i].size();
      for (int j = 0; j < polygons[i].size(); j++) {
        vertices.push_back(cdt.insert(pl.to_2d(m_lcc.point(m_lcc.template dart_of_attribute<0>(polygons[i][j])))));
        auto it = va2vh.insert(std::make_pair(polygons[i][j], vertices.size() - 1));
        CGAL_assertion(it.second);

        vertices.back()->info().i = i;
        vertices.back()->info().j = j;
        vertices.back()->info().p = pl.to_2d(m_lcc.point(m_lcc.template dart_of_attribute<0>(polygons[i][j])));
        vertices.back()->info().dh = polygons[i][j];

        if (j >= 1)
          cdt.insert_constraint(vertices[vertices.size() - 2], vertices.back());
      }
      cdt.insert_constraint(vertices.back(), vertices[vertices.size() - polygons[i].size()]);
    }
    // check infinitive edges for outer polygon
    int outer = -1;
    auto& e = *(cdt.incident_edges(cdt.infinite_vertex()));
    auto a = e.first->vertex((e.second + 1) % 3);
    auto b = e.first->vertex((e.second + 2) % 3);

    if (a == cdt.infinite_vertex())
      outer = b->info().i;
    else
      outer = a->info().i;

    CGAL_assertion(outer != -1);

    // Distance matrix
    std::vector<FT> dist(polygons.size()* polygons.size(), (std::numeric_limits<FT>::max)());
    std::vector<std::pair<std::size_t, std::size_t> > closest_pts(polygons.size() * polygons.size(), std::make_pair(-1, -1));

    for (auto& edge : cdt.finite_edges()) {
      auto v1 = edge.first->vertex((edge.second + 1) % 3);
      auto v2 = edge.first->vertex((edge.second + 2) % 3);

      if (v1->info().i != v2->info().i) {
        std::size_t idx;
        if (v1->info().i == -1 || v2->info().i == -1)
          continue;
        if (v1->info().i < v2->info().i)
          idx = v1->info().i * polygons.size() + v2->info().i;
        else
          idx = v2->info().i * polygons.size() + v1->info().i;

        FT d = from_exact((v1->info().p - v2->info().p).squared_length());
        if (dist[idx] > d) {
          dist[idx] = d;
          closest_pts[idx] = std::make_pair(v1->info().dh, v2->info().dh);
        }
      }
    }

    std::vector<bool> merged(polygons.size(), false);
    for (std::size_t i = 0; i < polygons.size(); i++) {
      if (i == outer)
        continue;

      std::size_t idx;
      if (i < outer)
        idx = i * polygons.size() + outer;
      else
        idx = outer * polygons.size() + i;

      // For now merge all polygons into outer if possible
      if (dist[idx] < (std::numeric_limits<FT>::max)()) {
        std::size_t in_target, in_source;
        for (in_target = 0; in_target < polygons[outer].size(); in_target++)
          if (polygons[outer][in_target] == closest_pts[idx].first || polygons[outer][in_target] == closest_pts[idx].second)
            break;

        for (in_source = 0; in_source < polygons[i].size(); in_source++)
            if (polygons[i][in_source] == closest_pts[idx].first || polygons[i][in_source] == closest_pts[idx].second)
              break;

        std::size_t former_end = polygons[outer].size() - 1;

        polygons[outer].resize(polygons[outer].size() + polygons[i].size() + 2);

        for (std::size_t j = 0; j != former_end - in_target + 1; j++)
          polygons[outer][polygons[outer].size() - j - 1] = polygons[outer][former_end - j];

        for (std::size_t j = 0; j < polygons[i].size() + 1; j++) {
          std::size_t idx = (in_source + j) % polygons[i].size();
          polygons[outer][in_target + j + 1] = polygons[i][idx];
        }
      }
      else {
        std::cout << "ghost edge could not be placed" << std::endl;
        // Do I need a minimum spanning tree? https://www.boost.org/doc/libs/1_75_0/libs/graph/example/kruskal-example.cpp
      }
      polygons[i].clear();
    }
    if (outer != 0)
      polygons[0] = std::move(polygons[outer]);
    polygons.resize(1);
  }

  typename LCC::Dart_descriptor circulate_vertex_2d(typename LCC::Dart_descriptor dh) {
    CGAL_assertion(!is_border_edge(dh));

    const Face_attribute& fa = m_lcc.template attribute<2>(dh);
    auto& finfo = m_lcc.template info_of_attribute<2>(fa);

    typename LCC::Dart_descriptor dh2 = m_lcc.template beta<2>(dh);

    std::size_t idx = 1;

    do {
      Face_attribute fa2 = m_lcc.template attribute<2>(dh2);
      auto& finfo2 = m_lcc.template info_of_attribute<2>(fa2);
      if (finfo2.input_polygon_index == finfo.input_polygon_index) {
        CGAL_assertion(fa != fa2);
        if (finfo2.input_polygon_index == -7) {
          if (finfo2.plane == finfo.plane || finfo2.plane == finfo.plane.opposite())
            return dh2;
        }
        else return dh2;
      }
      dh2 = m_lcc.template beta<3, 2>(dh2);
      idx++;

    } while (dh2 != dh);

    // dh is a border edge
    CGAL_assertion(false);

    return dh2;
  }

  void collect_border(typename LCC::Dart_descriptor dh, std::vector<bool>& processed, std::vector<std::vector<std::size_t> >& borders) {
    processed[dh] = true;

    if (!m_labels[m_lcc.template info<3>(dh).volume_id + 6] == 1)
      std::cout << "collect_border called on dart of outside volume, dh " << dh << " volume_id " << m_lcc.template info<3>(dh).volume_id << std::endl;

    std::vector<std::size_t> border;
    border.push_back(m_lcc.template attribute<0>(dh));

    const Face_attribute& fa = m_lcc.template attribute<2>(dh);
    auto& finfo = m_lcc.template info_of_attribute<2>(fa);

    typename LCC::Dart_descriptor cur = dh;
    cur = m_lcc.template beta<1>(cur);

    std::size_t idx = 0;

    do {
      if (is_border_edge(cur)) {
        CGAL_assertion(!processed[cur]);
        processed[cur] = true;
        border.push_back(m_lcc.template attribute<0>(cur));

        if (!m_labels[m_lcc.template info<3>(cur).volume_id + 6] == 1)
          std::cout << "border collected from dart of outside volume, dh " << cur << " volume_id " << m_lcc.template info<3>(cur).volume_id << std::endl;
      }
      else
        cur = circulate_vertex_2d(cur);
      cur = m_lcc.template beta<1>(cur);
      idx++;
    } while(cur != dh);

    borders.push_back(std::move(border));
  }

  void write_face(const typename LCC::Dart_descriptor dh, const std::string& fn) {
    std::vector<Point_3> face;
    From_exact from_exact;

    Dart_descriptor n = dh;
    do {
      face.push_back(from_exact(m_lcc.point(n)));
      n = m_lcc.beta(n, 1);
    } while (n != dh);

    KSP_3::internal::dump_polygon(face, fn);
  }

  void write_edge(typename LCC::Dart_descriptor dh, const std::string& fn) {
    From_exact from_exact;
    std::ofstream vout(fn);
    vout << "2 " << from_exact(m_lcc.point(dh)) << " " << from_exact(m_lcc.point(m_lcc.beta<1>(dh))) << std::endl;
    vout.close();
  }

  void write_border(std::vector<std::size_t> &border, const std::string& fn) {
    From_exact from_exact;
    std::ofstream vout(fn);
    vout << (border.size() + 1);
    for (std::size_t k = 0; k < border.size(); k++) {
      vout << " " << from_exact(m_lcc.point(m_lcc.dart_of_attribute<0>(border[k])));
    }
    vout << " " << from_exact(m_lcc.point(m_lcc.dart_of_attribute<0>(border[0]))) << std::endl;
    vout.close();
  }

  void collect_connected_border(std::vector<std::vector<std::size_t> >& borders, const std::vector<int> &region_index, std::vector<std::vector<std::size_t> > &borders_per_region) {
    // Start extraction of a border from each dart (each dart is a 1/n-edge)
    // Search starting darts by searching faces
    //borders contains Attribute<0> handles casted to std::size_t
    std::vector<bool> processed(m_lcc.number_of_darts(), false);

    borders_per_region.resize(region_index.size());

    for (std::size_t i = 0;i<region_index.size();i++) {
      if (region_index[i] == -1)
        continue;

      typename LCC::Dart_descriptor dh = m_faces_lcc[i];

      Volume_attribute va = m_lcc.template attribute<3>(dh);
      const Face_attribute &fa = m_lcc.template attribute<2>(dh);
      auto finfo = m_lcc.template info_of_attribute<2>(fa);
      const auto& n = m_face_neighbors_lcc[m_attrib2index_lcc[fa]];

      // Belongs to reconstruction?
      if (m_labels[n.first] == m_labels[n.second]) {
        std::cout << "face skipped" << std::endl;
        continue;
      }

      std::size_t num_edges = m_lcc.template one_dart_per_incident_cell<1, 2>(dh).size();

      typename LCC::Dart_descriptor dh2 = dh;

      do {
        if (va != m_lcc.template attribute<3>(dh2)) {
          std::cout << "volume attribute mismatch" << std::endl;
        }

        if (!processed[dh2] && is_border_edge(dh2)) {
          borders_per_region[region_index[fa]].push_back(borders.size());

          collect_border(dh2, processed, borders);
        }
        dh2 = m_lcc.template beta<1>(dh2);
      } while (dh2 != dh);
    }
  }

  void collect_points_for_faces_lcc() {
    FT total_area = 0;
    m_total_inliers = 0;
    From_exact from_exact;

    std::vector<std::vector<Dart_descriptor> > poly2faces(m_kinetic_partition.input_planes().size());
    std::vector<Dart_descriptor> other_faces;
    for (auto& d : m_lcc.template one_dart_per_cell<2>()) {
      Dart_descriptor dh = m_lcc.dart_descriptor(d);
      if (m_lcc.template info<2>(dh).input_polygon_index >= 0)
        poly2faces[m_lcc.template info<2>(dh).input_polygon_index].push_back(dh);
      else
        other_faces.push_back(dh); // Contains faces originating from the octree decomposition as well as bbox faces
    }

    assert(m_kinetic_partition.input_planes().size() == m_regions.size());

    std::size_t next = 0, step = 1;
    for (std::size_t i = 0; i < m_kinetic_partition.input_planes().size(); i++) {

      std::vector<std::pair<Dart_descriptor, std::vector<std::size_t> > > mapping;

      std::vector<Point_3> pts;
      pts.reserve(m_regions[i].second.size());

      for (std::size_t j = 0; j < m_regions[i].second.size(); j++)
        pts.emplace_back(get(m_point_map, m_regions[i].second[j]));

      map_points_to_faces(i, pts, mapping);

      // Remap from mapping to m_face_inliers
      for (auto p : mapping) {
        const Face_attribute& f = m_lcc.template attribute<2>(p.first);
        std::size_t id = m_attrib2index_lcc[f];
        assert(m_face_inliers[id].size() == 0);

        m_face_inliers[m_attrib2index_lcc[m_lcc.template attribute<2>(p.first)]].resize(p.second.size());
        for (std::size_t k = 0; k < p.second.size(); k++)
          m_face_inliers[m_attrib2index_lcc[m_lcc.template attribute<2>(p.first)]][k] = m_regions[i].second[p.second[k]];

        m_total_inliers += p.second.size();
      }

      Plane_3 pl = from_exact(m_kinetic_partition.input_planes()[i]);

      for (std::size_t j = 0; j < poly2faces[i].size(); j++) {
        std::size_t idx = m_attrib2index_lcc[m_lcc.template attribute<2>(poly2faces[i][j])];
        m_face_area_lcc[idx] = 0;

        //multiple regions per input polygon

        Delaunay_2 tri;

        Dart_descriptor n = poly2faces[i][j];
        do {
          tri.insert(pl.to_2d(from_exact(m_lcc.point(n))));
          n = m_lcc.beta(n, 0);
        } while (n != poly2faces[i][j]);

        // Get area
        for (auto fit = tri.finite_faces_begin(); fit != tri.finite_faces_end(); ++fit) {
          const Triangle_2 triangle(
            fit->vertex(0)->point(),
            fit->vertex(1)->point(),
            fit->vertex(2)->point());
          m_face_area_lcc[idx] += triangle.area();
        }

        total_area += m_face_area_lcc[idx];
      }
    }

    // Handling face generated by the octree partition. They are not associated with an input polygon.
    for (std::size_t i = 0; i < other_faces.size(); i++) {
      std::vector<Point_3> face;
      std::size_t idx = m_attrib2index_lcc[m_lcc.template attribute<2>(other_faces[i])];

      Dart_descriptor n = other_faces[i];
      do {
        face.push_back(from_exact(m_lcc.point(n)));
        n = m_lcc.beta(n, 0);
      } while (n != other_faces[i]);

      Plane_3 pl;
      CGAL::linear_least_squares_fitting_3(face.begin(), face.end(), pl, CGAL::Dimension_tag<0>());

      Delaunay_2 tri;
      for (const Point_3& p : face)
        tri.insert(pl.to_2d(p));

      // Get area
      for (auto fit = tri.finite_faces_begin(); fit != tri.finite_faces_end(); ++fit) {
        const Triangle_2 triangle(
          fit->vertex(0)->point(),
          fit->vertex(1)->point(),
          fit->vertex(2)->point());
        m_face_area_lcc[idx] += triangle.area();
      }

      total_area += m_face_area_lcc[idx];
    }

    m_face_area_lcc.resize(m_faces_lcc.size(), 0);

    for (std::size_t i = 0; i < m_faces_lcc.size(); i++)
      m_face_area_lcc[i] = m_face_area_lcc[i] * 2.0 * m_total_inliers / total_area;
  }

  void count_volume_votes_lcc() {
    const int debug_volume = -1;
    FT total_volume = 0;
    std::size_t num_volumes = m_kinetic_partition.number_of_volumes();
    m_volume_votes.clear();
    m_volume_votes.resize(num_volumes, std::make_pair(0, 0));

    m_volumes.resize(num_volumes, 0);

    for (std::size_t i = 6; i < num_volumes; i++) {
      m_cost_matrix[0][i] = m_cost_matrix[1][i] = 0;
      m_volumes[i] = 0;
    }

    std::size_t count_faces = 0;
    std::size_t count_points = 0;

    From_exact from_exact;

    std::size_t idx = 0;

    for (std::size_t i = 0; i < m_faces_lcc.size(); i++) {
      std::size_t v[] = { std::size_t(-1), std::size_t(-1) };
      Point_3 c[2];
      std::size_t in[] = {0, 0}, out[] = {0, 0};

      std::size_t idx = 0;
      for (auto& vd : m_lcc.template one_dart_per_incident_cell<3, 2>(m_faces_lcc[i])) {
        typename LCC::Dart_descriptor vdh = m_lcc.dart_descriptor(vd);
        v[idx] = m_lcc.template info<3>(vdh).volume_id;
        c[idx] = from_exact(m_lcc.template info<3>(vdh).barycenter);
        idx++;
      }

      for (std::size_t p : m_face_inliers[i]) {
        const auto& point = get(m_point_map, p);
        const auto& normal = get(m_normal_map, p);

        count_points++;

        for (std::size_t j = 0; j < idx; j++) {
          const Vector_3 vec(point, c[j]);
          const FT dot_product = vec * normal;
          if (dot_product < FT(0))
            in[j]++;
          else
            out[j]++;
        }
      }

      for (std::size_t j = 0; j < idx; j++) {
        m_volume_votes[v[j]].first += in[j];
        m_volume_votes[v[j]].second += out[j];
        m_cost_matrix[0][v[j] + 6] += static_cast<double>(in[j]);
        m_cost_matrix[1][v[j] + 6] += static_cast<double>(out[j]);
      }
    }

    for (auto &d : m_lcc.template one_dart_per_cell<3>()) {
      typename LCC::Dart_descriptor dh = m_lcc.dart_descriptor(d);

      std::vector<Point_3> volume_vertices;

      std::size_t volume_index = m_lcc.template info<3>(dh).volume_id;

      // Collect all vertices of volume to calculate volume
      for (auto &fd : m_lcc.template one_dart_per_incident_cell<2, 3>(dh)) {
        typename LCC::Dart_descriptor fdh = m_lcc.dart_descriptor(fd);

        for (const auto &vd : m_lcc.template one_dart_per_incident_cell<0, 2>(fdh))
          volume_vertices.push_back(from_exact(m_lcc.point(m_lcc.dart_descriptor(vd))));
      }

      Delaunay_3 tri;
      for (const Point_3& p : volume_vertices)
        tri.insert(p);

      m_volumes[volume_index] = FT(0);
      for (auto cit = tri.finite_cells_begin(); cit != tri.finite_cells_end(); ++cit) {
        const auto& tet = tri.tetrahedron(cit);
        m_volumes[volume_index] += tet.volume();
      }

      total_volume += m_volumes[volume_index];
    }

    // Normalize volumes
    for (FT& v : m_volumes)
      v /= total_volume;

//     for (std::size_t i = 0; i < m_volumes.size(); i++)
//       std::cout << i << ": " << m_cost_matrix[0][i] << " o: " << m_cost_matrix[1][i] << " v: " << m_volumes[i] << std::endl;
  }

  template<typename NamedParameters>
  void create_planar_shapes(const NamedParameters& np) {

    if (m_points.size() < 3) {
      if (m_verbose) std::cout << "* no points found, skipping" << std::endl;
      return;
    }
    if (m_verbose) std::cout << "* getting planar shapes using region growing" << std::endl;

    // Parameters.
    const std::size_t k = parameters::choose_parameter(parameters::get_parameter(np, internal_np::k_neighbors), 12);
    const FT max_distance_to_plane = parameters::choose_parameter(parameters::get_parameter(np, internal_np::maximum_distance), FT(1));
    const FT max_accepted_angle = parameters::choose_parameter(parameters::get_parameter(np, internal_np::maximum_angle), FT(15));
    const std::size_t min_region_size = parameters::choose_parameter(parameters::get_parameter(np, internal_np::minimum_region_size), 50);

    m_detection_distance_tolerance = max_distance_to_plane;

    // Region growing.
    Neighbor_query neighbor_query = CGAL::Shape_detection::Point_set::make_k_neighbor_query(
      m_points, CGAL::parameters::k_neighbors(k));

    Region_type region_type = CGAL::Shape_detection::Point_set::make_least_squares_plane_fit_region(
      m_points,
      CGAL::parameters::
      maximum_distance(max_distance_to_plane).
      maximum_angle(max_accepted_angle).
      minimum_region_size(min_region_size));

    Sorting sorting = CGAL::Shape_detection::Point_set::make_least_squares_plane_fit_sorting(m_points, neighbor_query);
    sorting.sort();

    Region_growing region_growing(
      m_points, sorting.ordered(), neighbor_query, region_type);
    region_growing.detect(std::back_inserter(m_regions));

    std::size_t unassigned = 0;
    region_growing.unassigned_items(m_points, boost::make_function_output_iterator([&](const auto&) { ++unassigned; }));

    std::vector<std::vector<Point_3> > polys_debug;

    for (std::size_t i = 0; i < m_regions.size(); i++) {

      Indices region;
      for (auto& j : m_regions[i].second)
        region.push_back(j);

      store_convex_hull_shape(region, m_regions[i].first, polys_debug);
      //KSR_3::dump_polygon(polys_debug[i], std::to_string(i) + "-detected-region.ply");
    }

    KSP_3::internal::dump_polygons(polys_debug, "detected-" + std::to_string(m_regions.size()) + "-polygons.ply");

    // Convert indices.
    m_planar_regions.clear();
    m_planar_regions.reserve(m_regions.size());

    // Copy planes for regularization.
    std::vector<Plane_3> planes(m_regions.size());
    for (std::size_t i = 0; i < m_regions.size(); i++)
      planes[i] = m_regions[i].first;

    auto range = m_regions | boost::adaptors::transformed([](typename Region_growing::Primitive_and_region& pr)->Plane_3& {return pr.first; });

    std::size_t num_shapes = m_regions.size();

    const bool regularize_axis_symmetry = parameters::choose_parameter(parameters::get_parameter(np, internal_np::regularize_axis_symmetry), false);
    const bool regularize_coplanarity = parameters::choose_parameter(parameters::get_parameter(np, internal_np::regularize_coplanarity), false);
    const bool regularize_orthogonality = parameters::choose_parameter(parameters::get_parameter(np, internal_np::regularize_orthogonality), false);
    const bool regularize_parallelism = parameters::choose_parameter(parameters::get_parameter(np, internal_np::regularize_parallelism), false);
    const FT angle_tolerance = parameters::choose_parameter(parameters::get_parameter(np, internal_np::angle_tolerance), 25);
    const FT maximum_offset = parameters::choose_parameter(parameters::get_parameter(np, internal_np::maximum_offset), 0.01);

    // Regularize detected planes.

    CGAL::Shape_regularization::Planes::regularize_planes(range, m_points,
      CGAL::parameters::plane_index_map(region_growing.region_map())
      .point_map(m_point_map)
      .regularize_axis_symmetry(regularize_axis_symmetry)
      .regularize_orthogonality(regularize_orthogonality)
      .regularize_parallelism(regularize_parallelism)
      .regularize_coplanarity(regularize_coplanarity)
      .maximum_angle(angle_tolerance)
      .maximum_offset(maximum_offset));

    polys_debug.clear();

    for (std::size_t i = 0; i < m_regions.size(); i++) {

      Indices region;
      for (auto& j : m_regions[i].second)
        region.push_back(j);

      store_convex_hull_shape(region, m_regions[i].first, polys_debug);
      //KSR_3::dump_polygon(polys_debug[i], std::to_string(i) + "-detected-region.ply");
    }

    KSP_3::internal::dump_polygons(polys_debug, "regularized-" + std::to_string(m_regions.size()) + "-polygons.ply");

    // Merge coplanar regions
    for (std::size_t i = 0; i < m_regions.size() - 1; i++) {
      for (std::size_t j = i + 1; j < m_regions.size(); j++) {
        if (m_regions[i].first == m_regions[j].first || m_regions[i].first.opposite() == m_regions[j].first) {
          std::move(m_regions[j].second.begin(), m_regions[j].second.end(), std::back_inserter(m_regions[i].second));
          m_regions.erase(m_regions.begin() + j);
          j--;
        }
      }
    }

    // Estimate ground plane by finding a low mostly horizontal plane
    std::vector<std::size_t> candidates;
    FT low_z_peak = (std::numeric_limits<FT>::max)();
    FT bbox_min[] = { (std::numeric_limits<FT>::max)(), (std::numeric_limits<FT>::max)(), (std::numeric_limits<FT>::max)() };
    FT bbox_max[] = { -(std::numeric_limits<FT>::max)(), -(std::numeric_limits<FT>::max)(), -(std::numeric_limits<FT>::max)() };
    for (const auto& p : m_points) {
      const auto& point = get(m_point_map, p);
      for (int i = 0; i < 3; i++) {
        bbox_min[i] = (std::min)(point[i], bbox_min[i]);
        bbox_max[i] = (std::max)(point[i], bbox_max[i]);
      }
    }

    FT bbox_center[] = { 0.5 * (bbox_min[0] + bbox_max[0]), 0.5 * (bbox_min[1] + bbox_max[1]), 0.5 * (bbox_min[2] + bbox_max[2]) };

    for (std::size_t i = 0; i < m_regions.size(); i++) {
      Vector_3 d = m_regions[i].first.orthogonal_vector();
      if (abs(d.z()) > 0.98) {
        candidates.push_back(i);
        FT z = m_regions[i].first.projection(Point_3(bbox_center[0], bbox_center[1], bbox_center[2])).z();
        low_z_peak = (std::min<FT>)(z, low_z_peak);
      }
    }

    m_ground_polygon_index = -1;
    std::vector<std::size_t> other_ground;
    for (std::size_t i = 0; i < candidates.size(); i++) {
      Vector_3 d = m_regions[candidates[i]].first.orthogonal_vector();
      FT z = m_regions[candidates[i]].first.projection(Point_3(bbox_center[0], bbox_center[1], bbox_center[2])).z();
      if (z - low_z_peak < max_distance_to_plane) {
        if (m_ground_polygon_index == -1)
          m_ground_polygon_index = candidates[i];
        else
          other_ground.push_back(candidates[i]);
      }
    }

    if (m_ground_polygon_index != -1) {

      for (std::size_t i = 0; i < other_ground.size(); i++)
        std::move(m_regions[other_ground[i]].second.begin(), m_regions[other_ground[i]].second.end(), std::back_inserter(m_regions[m_ground_polygon_index].second));

      if (m_verbose)
        std::cout << "ground polygon " << m_ground_polygon_index << ", merging other polygons";

      while (other_ground.size() != 0) {
        m_regions.erase(m_regions.begin() + other_ground.back());
        std::cout << " " << other_ground.back();
        other_ground.pop_back();
      }
      std::cout << std::endl;

      std::vector<Point_3> ground_plane;
      ground_plane.reserve(m_regions[m_ground_polygon_index].second.size());
      for (std::size_t i = 0; i < m_regions[m_ground_polygon_index].second.size(); i++) {
        ground_plane.push_back(get(m_point_map, m_regions[m_ground_polygon_index].second[i]));
      }

      CGAL::linear_least_squares_fitting_3(ground_plane.begin(), ground_plane.end(), m_regions[m_ground_polygon_index].first, CGAL::Dimension_tag<0>());

      if (m_regions[m_ground_polygon_index].first.orthogonal_vector().z() < 0)
        m_regions[m_ground_polygon_index].first = m_regions[m_ground_polygon_index].first.opposite();
    }

    polys_debug.clear();

    for (std::size_t i = 0; i < m_regions.size(); i++) {

      Indices region;
      for (auto& j : m_regions[i].second)
        region.push_back(j);

      store_convex_hull_shape(region, m_regions[i].first, polys_debug);
      //KSP_3::dump_polygon(polys_debug[i], std::to_string(i) + "-detected-region.ply");
    }

    if (m_debug)
      KSP_3::internal::dump_polygons(polys_debug, "merged-" + std::to_string(m_regions.size()) + "-polygons.ply");

    std::vector<Plane_3> pl;

    std::size_t idx = 0;
    for (const auto& p : m_regions) {
      bool exists = false;
      for (std::size_t i = 0; i < pl.size(); i++)
        if (pl[i] == p.first || pl[i].opposite() == p.first)
          exists = true;

      if (!exists)
        pl.push_back(p.first);

      idx++;
    }

    for (const auto& pair : m_regions) {
      Indices region;
      for (auto& i : pair.second)
        region.push_back(i);
      m_planar_regions.push_back(region);

      const std::size_t shape_idx = add_convex_hull_shape(region, pair.first);
      CGAL_assertion(shape_idx != std::size_t(-1));
      m_region_map[shape_idx] = region;
    }
    CGAL_assertion(m_planar_regions.size() == m_regions.size());

    if (m_verbose) {
      std::cout << "found " << num_shapes << " planar shapes regularized into " << m_planar_regions.size() << std::endl;
      std::cout << "from " << m_points.size() << " input points " << unassigned << " remain unassigned" << std::endl;
    }
  }

  void map_points_to_faces(const std::size_t polygon_index, const std::vector<Point_3>& pts, std::vector<std::pair<typename LCC::Dart_descriptor, std::vector<std::size_t> > >& face_to_points) {
    std::vector<Index> faces;

    if (polygon_index >= m_kinetic_partition.input_planes().size())
      assert(false);

    From_exact from_exact;

    const typename Intersection_kernel::Plane_3& pl = m_kinetic_partition.input_planes()[polygon_index];
    const Plane_3 inexact_pl = from_exact(pl);
    std::vector<Point_2> pts2d;
    pts2d.reserve(pts.size());

    for (const Point_3& p : pts)
      pts2d.push_back(inexact_pl.to_2d(p));

    // Iterate over all faces of the lcc
    for (Dart& d : m_lcc.template one_dart_per_cell<2>()) {
      Dart_descriptor dd = m_lcc.dart_descriptor(d);
      if (m_lcc.template info<2>(m_lcc.dart_descriptor(d)).input_polygon_index != polygon_index || !m_lcc.template info<2>(m_lcc.dart_descriptor(d)).part_of_initial_polygon)
        continue;

      // No filtering of points per partition

      face_to_points.push_back(std::make_pair(m_lcc.dart_descriptor(d), std::vector<std::size_t>()));

      auto& info = m_lcc.template info<2>(m_lcc.dart_descriptor(d));

      std::vector<Point_2> vts2d;
      vts2d.reserve(m_lcc.template one_dart_per_incident_cell<0, 2>(m_lcc.dart_descriptor(d)).size());

      typename LCC::Dart_descriptor n = dd;
      do {
        vts2d.push_back(inexact_pl.to_2d(from_exact(m_lcc.point(n))));
        n = m_lcc.beta(n, 0);
      } while (n != dd);

      Polygon_2<Kernel> poly(vts2d.begin(), vts2d.end());
      if (poly.is_clockwise_oriented())
        std::reverse(vts2d.begin(), vts2d.end());

      for (std::size_t i = 0; i < pts2d.size(); i++) {
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

        face_to_points.back().second.push_back(i);
      }
    }
  }

/*
  const Plane_3 fit_plane(const std::vector<std::size_t>& region) const {

    std::vector<Point_3> points;
    points.reserve(region.size());
    for (const std::size_t idx : region) {
      CGAL_assertion(idx < m_points.size());
      points.push_back(get(m_point_map, idx));
    }
    CGAL_assertion(points.size() == region.size());

    Plane_3 fitted_plane;
    Point_3 fitted_centroid;
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
*/

  void set_outside_volumes(bool ground, std::vector<std::vector<double> >& cost_matrix) const {
    // Setting preferred outside label for bbox plane nodes
    // Order:
    // 0 zmin
    // 1 ymin
    // 2 xmax
    // 3 ymax
    // 4 xmin
    // 5 zmax
    const double force = static_cast<double>(m_total_inliers * 3);
    // 0 - cost for labelled as outside
    cost_matrix[0][0] = 0;
    cost_matrix[0][1] = 0;
    cost_matrix[0][2] = 0;
    cost_matrix[0][3] = 0;
    cost_matrix[0][4] = 0;
    cost_matrix[0][5] = 0;
    // 1 - cost for labelled as inside
    cost_matrix[1][0] = 0;
    cost_matrix[1][1] = 0;
    cost_matrix[1][2] = 0;
    cost_matrix[1][3] = 0;
    cost_matrix[1][4] = 0;
    cost_matrix[1][5] = 0;

    if (m_ground_polygon_index != static_cast<std::size_t>(-1) && ground) {
      if (m_verbose)
        std::cout << "using estimated ground plane for reconstruction" << std::endl;
      cost_matrix[0][0] = force;
      cost_matrix[0][1] = 0;
      cost_matrix[0][2] = 0;
      cost_matrix[0][3] = 0;
      cost_matrix[0][4] = 0;
      cost_matrix[1][0] = 0;
      cost_matrix[1][1] = force;
      cost_matrix[1][2] = force;
      cost_matrix[1][3] = force;
      cost_matrix[1][4] = force;
    }
  }
};

} // namespace CGAL


#endif // CGAL_KINETIC_SURFACE_RECONSTRUCTION_3_H
