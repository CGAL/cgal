// Copyright (c) 2019 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Simon Giraudot, Dmitry Anisimov, Sven Oesau

#ifndef CGAL_KINETIC_SHAPE_RECONSTRUCTION_3_H
#define CGAL_KINETIC_SHAPE_RECONSTRUCTION_3_H

//#include <CGAL/license/Kinetic_shape_reconstruction.h>
#include <CGAL/Kinetic_shape_partitioning_3.h>
#include <CGAL/KSR/debug.h>

namespace CGAL
{
/*!
* \ingroup PkgKineticPartition
  \brief Piece-wise linear reconstruction via inside/outside labeling of a kinetic partition using graph cut.

  \tparam GeomTraits
    must be a model of `KineticShapePartitionTraits_3`.

  \tparam NormalMap
    must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`. It must map the elements in `KineticShapePartitionTraits_3::Input_range` to `Vector_3`.
*/
template<typename Traits, typename NormalMap>
class Kinetic_shape_reconstruction_3 {
public:
  using Input_range = typename Traits::Input_range;
  using Kernel = typename Traits::Kernel;
  using Intersection_Kernel = typename Traits::Intersection_Kernel;

  using FT = typename Kernel::FT;

  using Point_2 = typename Traits::Point_2;
  using Point_3 = typename Traits::Point_3;
  using Plane_3 = typename Traits::Plane_3;

  using Point_map = typename Traits::Point_map;
  using Normal_map = NormalMap;

  using Indices = std::vector<std::size_t>;
  using Polygon_3 = std::vector<Point_3>;

  using KSP = Kinetic_shape_partitioning_3<Traits>;

  using Mesh = Surface_mesh<Point_3>;

  using Neighbor_query_3 = CGAL::Shape_detection::Point_set::K_neighbor_query<Kernel, Input_range, Point_map>;
  using Planar_region = CGAL::Shape_detection::Point_set::Least_squares_plane_fit_region<Kernel, Input_range, Point_map, Normal_map>;
  using Planar_sorting = CGAL::Shape_detection::Point_set::Least_squares_plane_fit_sorting<Kernel, Input_range, Neighbor_query_3, Point_map>;
  using Region_growing_3 = CGAL::Shape_detection::Region_growing<Input_range, Neighbor_query_3, Planar_region, typename Planar_sorting::Seed_map>;
  /*!
  \brief Creates a Kinetic_shape_reconstruction_3 object.

  \param verbose
  provides information on std::cout. The default is false.

  \param debug
  writes intermediate results into ply files. The default is false.

  */
  Kinetic_shape_reconstruction_3(const Input_range &input_range, bool verbose = false, bool debug = false) : m_kinetic_partition(verbose, debug), m_points(input_range) {}

  /*!
    \brief Detects shapes in the provided point cloud

    \tparam InputRange
    must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam NamedParameters
    a sequence of \ref bgl_namedparameters "Named Parameters"

    \param input_range
    an instance of `InputRange` with 3D points and corresponding 3D normal vectors

    \param np
    an instance of `NamedParameters`.

  \cgalNamedParamsBegin
    \cgalParamNBegin{point_map}
      \cgalParamDescription{a property map associating points to the elements of `input_range`}
      \cgalParamDefault{`PointMap()`}
    \cgalParamNEnd
    \cgalParamNBegin{normal_map}
      \cgalParamDescription{a property map associating normals to the elements of `input_range`}
      \cgalParamDefault{`NormalMap()`}
    \cgalParamNEnd
    \cgalParamNBegin{k_neighbors}
      \cgalParamDescription{the number of returned neighbors per each query point}
      \cgalParamType{`std::size_t`}
      \cgalParamDefault{12}
    \cgalParamNEnd
    \cgalParamNBegin{distance_threshold}
      \cgalParamDescription{maximum distance from a point to a planar shape}
      \cgalParamType{`GeomTraits::FT`}
      \cgalParamDefault{1}
    \cgalParamNEnd
    \cgalParamNBegin{angle_threshold}
      \cgalParamDescription{maximum angle in degrees between the normal of a point and the plane normal}
      \cgalParamType{`GeomTraits::FT`}
      \cgalParamDefault{25 degrees}
    \cgalParamNEnd
    \cgalParamNBegin{min_region_size}
      \cgalParamDescription{minimum number of 3D points a region must have}
      \cgalParamType{`std::size_t`}
      \cgalParamDefault{5}
    \cgalParamNEnd
  \cgalNamedParamsEnd

  */
  template<
    typename InputRange,
    typename CGAL_NP_TEMPLATE_PARAMETERS>
  std::size_t detect_planar_shapes(
    InputRange input_range,
    const CGAL_NP_CLASS& np = parameters::default_values()) {

    if (m_verbose)
      std::cout << std::endl << "--- DETECTING PLANAR SHAPES: " << std::endl;

    m_planes.clear();
    m_polygons.clear();
    m_region_map.clear();

    m_point_map = Point_set_processing_3_np_helper<InputRange, CGAL_NP_CLASS, Point_map>::get_point_map(input_range, np);
    m_normal_map = Point_set_processing_3_np_helper<InputRange, CGAL_NP_CLASS, Normal_map>::get_normal_map(input_range, np);

    create_planar_shapes(np);

    CGAL_assertion(m_planes.size() == m_polygons.size());
    CGAL_assertion(m_polygons.size() == m_region_map.size());

    return m_polygons.size();
  }

  /*!
    \brief Regularizes detected planar shapes by using `CGAL::Shape_regularization::Planes::regularize_planes` and merging coplanar planes afterwards.

    \tparam NamedParameters
    a sequence of \ref bgl_namedparameters "Named Parameters"

    \param np
    an instance of `NamedParameters`.

    \cgalNamedParamsBegin
    \cgalParamNBegin{point_map}
      \cgalParamDescription{a property map associating points to the elements of `input_range` that has been passed to detect_planar_shapes}
      \cgalParamType{bool}
      \cgalParamDefault{false}
    \cgalParamNEnd
    \cgalParamNBegin{normal_map}
      \cgalParamDescription{a property map associating normals to the elements of `input_range` that has been passed to detect_planar_shapes}
      \cgalParamType{bool}
      \cgalParamDefault{false}
    \cgalParamNEnd
    \cgalParamNBegin{maximum_angle}
      \cgalParamDescription{maximum allowed angle in degrees between plane normals used for parallelism, orthogonality, and axis symmetry}
      \cgalParamType{FT}
      \cgalParamDefault{25 degrees}
    \cgalParamNEnd
    \cgalParamNBegin{maximum_offset}
      \cgalParamDescription{maximum allowed orthogonal distance between two parallel planes such that they are considered to be coplanar}
      \cgalParamType{FT}
      \cgalParamDefault{0.01}
    \cgalParamNEnd
    \cgalParamNBegin{regularize_parallelism}
      \cgalParamDescription{indicates whether parallelism should be regularized or not}
      \cgalParamType{bool}
      \cgalParamDefault{true}
    \cgalParamNEnd
    \cgalParamNBegin{regularize_orthogonality}
      \cgalParamDescription{indicates whether orthogonality should be regularized or not}
      \cgalParamType{bool}
      \cgalParamDefault{true}
    \cgalParamNEnd
    \cgalParamNBegin{regularize_coplanarity}
      \cgalParamDescription{indicates whether coplanarity should be regularized or not}
      \cgalParamType{bool}
      \cgalParamDefault{true}
    \cgalParamNEnd
  \cgalNamedParamsEnd

  */
  template<typename NamedParameters>
  std::size_t regularize_shapes(
    const NamedParameters& np) {

    /*if (m_verbose)
      std::cout << std::endl << "--- REGULARIZING PLANAR SHAPES: " << std::endl;

    const bool regularize = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::regularize), false);
    if (!regularize) {
      if (m_verbose) std::cout << "* user-defined, skipping" << std::endl;
      return true;
    }

    if (m_polygons.size() == 0) {
      if (m_verbose) std::cout << "* no planes found, skipping" << std::endl;
      return false;
    }

    // Regularize.

    std::vector<Plane_3> planes;
    std::vector<Indices> regions;
    create_planes_and_regions(planes, regions);

    CGAL_assertion(planes.size() > 0);
    CGAL_assertion(planes.size() == regions.size());

    Plane_map plane_map;
    Point_to_plane_map point_to_plane_map(m_input_range, regions);
    CGAL::Shape_regularization::Planes::regularize_planes(
      m_input_range,
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
      const auto& plane = planes[i];
      const auto& region = regions[i];

      const std::size_t shape_idx = add_planar_shape(region, plane);
      CGAL_assertion(shape_idx != std::size_t(-1));
      m_region_map[shape_idx] = region;
    }
    CGAL_assertion(m_polygons.size() == num_polygons);
    CGAL_assertion(m_polygons.size() == m_planes.size());
    CGAL_assertion(m_polygons.size() == m_region_map.size());

    if (m_verbose)
      std::cout << "* num regularized planes: " << m_planes.size() << std::endl;


    if (m_debug)
      dump_polygons("regularized-planar-shapes");*/
    return true;
  }

  /*!
  \brief Retrieves the detected shapes.

  @return
  vector with a plane equation for each detected planar shape.

  \pre `successful shape detection`
  */
  const std::vector<Plane_3>& detected_planar_shapes() {
    return m_planes;
  }

  /*!
  \brief Retrieves the indices of detected shapes.

    @return
    indices into `input_range` for each detected planar shape in vectors.

    \pre `successful shape detection`
  */
  const std::vector<std::vector<std::size_t> >& detected_planar_shape_indices() {
    return m_planar_regions;
  }

  /*!
  \brief initializes the kinetic partitioning.

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
      \cgalParamDescription{Factor for extension of the bounding box of the input data to be used for the partitioning.}
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
      \cgalParamDefault{5}
    \cgalParamNEnd
  \cgalNamedParamsEnd

    \pre `successful shape detection`
  */
  template<typename CGAL_NP_TEMPLATE_PARAMETERS>
  bool initialize_partitioning(const CGAL_NP_CLASS& np = parameters::default_values()) {
    if (m_polygons.size() == 0) {
      std::cout << "No planar shapes available to create kinetic partitioning." << std::endl;
      return false;
    }

    using Polygon_3 = std::vector<Point_3>;
    using Polygon_map = CGAL::Identity_property_map<Polygon_3>;

    return m_kinetic_partition.initialize(m_polygons, Polygon_map(), np);
  }

  /*!
  \brief Propagates the kinetic polygons in the initialized partition.

  \param k
  maximum number of allowed intersections for each input polygon before its expansion stops.

  @return
  success of kinetic partitioning.

  \pre `successful initialization`
  */
  bool partition(std::size_t k) {
    return m_kinetic_partition.partition(k);
  }

  /*!
  \brief Access to the kinetic partitioning.

  @return
  created kinetic partitioning data structure

  \pre `successful partitioning`
  */
  const Kinetic_shape_partitioning_3<Traits>& partitioning() const {
    return m_kinetic_partition;
  }

  /*!
  \brief Creates the visibility (data-) and regularity energy terms from the input point cloud and the kinetic partitioning.

  @return
  success.

  \pre `successful initialization`
  */
  bool setup_energyterms() {
    if (m_kinetic_partition.number_of_volumes() == 0) {
      if (m_verbose) std::cout << "Kinetic partition is not constructed or does not have volumes" << std::endl;
      return false;
    }

/*
    if (m_verbose) std::cout << "* computing visibility ... ";
    std::map<std::size_t, Indices> face2points;
    assign_points_to_pfaces(face2points);
    const Visibility visibility(
      m_data, face2points, m_point_map_3, m_normal_map_3);

    CGAL_assertion(m_data.volumes().size() > 0);
    visibility.compute(m_data.volumes());
    //dump_visibility("visibility/visibility", pface_points);

    if (m_verbose) {
      std::cout << "done" << std::endl;
      std::cout << "* applying graphcut ... ";
    }

    const FT beta = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::graphcut_beta), FT(1) / FT(2));*/

    return false;
  }

  /*!
  \brief Provides the data and regularity energy terms for reconstruction via graph-cut.

  \param edges
  contains a vector of pairs of volume indices. Indicates which volumes should be connected in the graph cut formulation.

  \param edge_costs
  contains the cost for each edge specified in `edges` for two labels with different labels. For equal labels, the cost is 0. Needs to be index compatible to the `edges` parameter.

  \param cost_matrix
  provides the cost of a label for each volume cell. The first index corresponds to the label and the second index corresponds to the volume index.

  @return
  fails if the dimensions of parameters does not match the kinetic partitioning.

  \pre `successful initialization`
  */
  template<typename NamedParameters>
  bool setup_energyterms(
    const std::vector< std::pair<std::size_t, std::size_t> >& edges,
    const std::vector<double>& edge_costs,
    const std::vector< std::vector<double> >& cost_matrix);

  /*!
  \brief Uses graph-cut to solve an solid/empty labeling of the volumes of the kinetic partition.

  \param lambda
  trades the impact of the data term for impact of the regularization term. Should be in the range [0, 1).

  @return
  success of reconstruction.

  \pre `successful initialization`
  */
  bool reconstruct(FT lambda) {
    return false;
  }

  /*!
  \brief Provides the reconstructed surface mesh

  \param mesh
  a mesh object to store the reconstructed surface.

  \pre `successful reconstruction`
  */
  void output_reconstructed_model(Mesh& mesh);

private:
  bool m_verbose;
  bool m_debug;

  const Input_range &m_points;
  Point_map m_point_map;
  Normal_map m_normal_map;

  std::vector<std::vector<std::size_t> > m_planar_regions;
  std::map<std::size_t, Indices> m_region_map;

  std::vector<Plane_3> m_planes;
  std::vector<Polygon_3> m_polygons;
  KSP m_kinetic_partition;

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
    return shape_idx;
  }

  template<typename NamedParameters>
  void create_planar_shapes(const NamedParameters& np) {

    if (m_points.size() < 3) {
      if (m_verbose) std::cout << "* no points found, skipping" << std::endl;
      return;
    }
    if (m_verbose) std::cout << "* getting planar shapes using region growing" << std::endl;

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
    Neighbor_query_3 neighbor_query(m_points, k, m_point_map);

    Planar_region planar_region(m_points,
      max_distance_to_plane, max_accepted_angle, min_region_size,
      m_point_map, m_normal_map);

    Planar_sorting sorting(
      m_points, neighbor_query, m_point_map);
    sorting.sort();

    std::vector<Indices> result;
    Region_growing_3 region_growing(
      m_points, neighbor_query, planar_region, sorting.seed_map());
    region_growing.detect(std::back_inserter(result));

    // Convert indices.
    m_planar_regions.clear();
    m_planar_regions.reserve(result.size());

    Indices region;
    for (const auto& indices : result) {
      region.clear();
      for (const std::size_t index : indices) {
        region.push_back(index);
      }
      m_planar_regions.push_back(region);
      const auto plane = fit_plane(region);
      const std::size_t shape_idx = add_convex_hull_shape(region, plane);
      CGAL_assertion(shape_idx != std::size_t(-1));
      m_region_map[shape_idx] = region;
    }
    CGAL_assertion(m_planar_regions.size() == result.size());

    if (m_verbose)
      std::cout << "* found " << m_polygons.size() << " planar shapes" << std::endl;
  }

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
};


} // namespace CGAL


#endif // CGAL_KINETIC_SHAPE_RECONSTRUCTION_3_H
