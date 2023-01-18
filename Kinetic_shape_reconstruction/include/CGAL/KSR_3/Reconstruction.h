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

#ifndef CGAL_KSR_3_RECONSTRUCTION_H
#define CGAL_KSR_3_RECONSTRUCTION_H

// #include <CGAL/license/Kinetic_shape_reconstruction.h>

// CGAL includes.
#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_point_set.h>
#include <CGAL/Shape_regularization/regularize_planes.h>

#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

// Internal includes.
#include <CGAL/KSR/enum.h>
#include <CGAL/KSR/utils.h>
#include <CGAL/KSR/debug.h>
#include <CGAL/KSR/property_map.h>
#include <CGAL/KSR_3/Data_structure.h>
#include <CGAL/KSR_3/Visibility.h>
#include <CGAL/KSR_3/Graphcut.h>

namespace CGAL {
namespace KSR_3 {

#ifdef DOXYGEN_RUNNING
/*!
  \brief Piece-wise linear reconstruction via inside/outside labeling of a kinetic partition using graph cut.

  \tparam Kernel
    must be a model of `Kernel`.

  \tparam IntersectionKernel
    must be a model of `Kernel`. Is used for the creation of the intersection graph. An exact kernel is suggested.
*/
template<Kernel, Intersection_Kernel>
class Kinetic_reconstruction_3 {
  /*!
    \brief Creates the kinetic partitioning of the bounding box.

    \tparam InputRange
    must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam SemanticMap
    must be an `LvaluePropertyMap` whose key type is the value type of the input
    range and value type is `std::size_t`.

    \tparam NamedParameters
    a sequence of \ref bgl_namedparameters "Named Parameters"

    \param input_range
    an instance of `InputRange` with 3D points and corresponding 3D normal vectors

  \cgalNamedParamsBegin
    \cgalParamNBegin{point_map}
      \cgalParamDescription{}
      \cgalParamType{bool}
      \cgalParamDefault{false}
    \cgalParamNEnd
    \cgalParamNBegin{normal_map}
      \cgalParamDescription{}
      \cgalParamType{bool}
      \cgalParamDefault{false}
    \cgalParamNEnd
    \cgalParamNBegin{semantic_map}
      \cgalParamDescription{A `LvaluePropertyMap` whose key type is the value type of the input range and value type is `std::size_t`.}
      \cgalParamType{bool}
      \cgalParamDefault{false}
    \cgalParamNEnd
  \cgalNamedParamsEnd

  */
  template<
    typename InputRange,
    typename NamedParameters>
  std::size_t detect_planar_shapes(
    InputRange input_range,
    const NamedParameters& np);

  /*!
    \brief Regularizes detected planar shapes by using `CGAL::Shape_regularization::Planes::regularize_planes` and merging coplanar planes afterwards.

    \tparam NamedParameters
    a sequence of \ref bgl_namedparameters "Named Parameters"

    \cgalNamedParamsBegin
    \cgalParamNBegin{maximum_angle}
      \cgalParamDescription{}
      \cgalParamType{FT}
      \cgalParamDefault{25 degrees}
    \cgalParamNEnd
    \cgalParamNBegin{maximum_offset}
      \cgalParamDescription{}
      \cgalParamType{FT}
      \cgalParamDefault{0.01}
    \cgalParamNEnd
    \cgalParamNBegin{regularize_parallelism}
      \cgalParamDescription{}
      \cgalParamType{bool}
      \cgalParamDefault{true}
    \cgalParamNEnd
    \cgalParamNBegin{regularize_orthogonality}
      \cgalParamDescription{}
      \cgalParamType{bool}
      \cgalParamDefault{true}
    \cgalParamNEnd
    \cgalParamNBegin{regularize_coplanarity}
      \cgalParamDescription{}
      \cgalParamType{bool}
      \cgalParamDefault{true}
    \cgalParamNEnd
  \cgalNamedParamsEnd

  */
  template<typename NamedParameters>
  std::size_t regularize_shapes(
    const NamedParameters& np);

  /*!
  \brief Retrieves the detected shapes.

    \param indices
    will be used to store the indices into the input range for each detected planar shape.

    \param planes
    will be used to store the plane equation of each detected planar shape.

    \pre `successful shape detection`
  */
  void detected_shapes(std::vector<std::vector<std::size_t> >& indices, std::vector<Plane_3>& planes);

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

    \pre `successful shape detection`
  */
  template<typename NamedParameters>
  bool initialize_partitioning(const NamedParameters& np);

  /*!
  \brief Propagates the kinetic polygons in the initialized partition.

  \param k
  maximum number of allowed intersections for each input polygon before its expansion stops.

  @return
  success of kinetic partitioning.

  \pre `successful initialization`
  */
  bool partitioning(std::size_t k);

  /*!
  \brief Access to the kinetic partitioning.

  @return
  created kinetic partitioning data structure

  \pre `successful partitioning`
  */
  const Kinetic_partitioning_3<Kernel, Intersection_Kernel>& get_partitioning() const;

  /*!
  \brief Creates the visibility (data-) and regularity energy terms from the input point cloud and the kinetic partitioning.

  @return
  success.

  \pre `successful initialization`
  */
  template<typename NamedParameters>
  bool setup_energyterms(const NamedParameters& np);

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
    const std::vector< std::vector<double> >& cost_matrix,
    const NamedParameters& np);

  /*!
  \brief Propagates the kinetic polygons in the initialized partition.

  \param lambda
  trades the impact of the data term for impact of the regularization term. Should be in the range [0, 1).

  @return
  success of reconstruction.

  \pre `successful initialization`
  */
  bool reconstruct(FT lambda);

  /*!
  \brief Propagates the kinetic polygons in the initialized partition.

  \param lambda
  trades the impact of the data term for impact of the regularization term. Should be in the range [0, 1).

  @return
  success of reconstruction.

  \pre `successful initialization`
  */
  void output_reconstructed_model(Polygon_mesh& polygon_mesh);
}
#else

template<
typename InputRange,
typename PointMap,
typename VectorMap,
typename SemanticMap,
typename GeomTraits>
class Reconstruction {

public:
  using Input_range  = InputRange;
  using Point_map    = PointMap;
  using Vector_map   = VectorMap;
  using Semantic_map = SemanticMap;
  using Kernel       = GeomTraits;

private:
  using FT        = typename Kernel::FT;
  using Point_2   = typename Kernel::Point_2;
  using Line_2    = typename Kernel::Line_2;
  using Point_3   = typename Kernel::Point_3;
  using Plane_3   = typename Kernel::Plane_3;
  using Vector_2  = typename Kernel::Vector_2;
  using Vector_3  = typename Kernel::Vector_3;
  using Segment_2 = typename Kernel::Segment_2;
  using Segment_3 = typename Kernel::Segment_3;

  using Data_structure = KSR_3::Data_structure<Kernel>;
  using PFace          = typename Data_structure::PFace;

  using Point_map_3        = KSR::Item_property_map<Input_range, Point_map>;
  using Vector_map_3       = KSR::Item_property_map<Input_range, Vector_map>;
  using Plane_map          = CGAL::Identity_property_map<Plane_3>;
  using Point_to_plane_map = CGAL::Shape_detection::RG::Point_to_shape_index_map;

  using Semantic_label    = KSR::Semantic_label;
  using Planar_shape_type = KSR::Planar_shape_type;

  using Indices     = std::vector<std::size_t>;
  using Polygon_3   = std::vector<Point_3>;
  using Polygon_map = CGAL::Identity_property_map<Polygon_3>;

  using IK       = CGAL::Exact_predicates_inexact_constructions_kernel;
  using IPoint_2 = typename IK::Point_2;
  using ILine_2  = typename IK::Line_2;
  using IPoint_3 = typename IK::Point_3;
  using IPlane_3 = typename IK::Plane_3;

  using Converter = CGAL::Cartesian_converter<Kernel, IK>;
  using From_EK = CGAL::Cartesian_converter<CGAL::Exact_predicates_exact_constructions_kernel, IK>;

  struct Vertex_info { FT z = FT(0); };
  struct Face_info { };

  using Fbi = CGAL::Triangulation_face_base_with_info_2<Face_info, Kernel>;
  using Fb  = CGAL::Alpha_shape_face_base_2<Kernel, Fbi>;

  using Vbi = CGAL::Triangulation_vertex_base_with_info_2<Vertex_info, Kernel>;
  using Vb  = CGAL::Alpha_shape_vertex_base_2<Kernel, Vbi>;

  using Tds         = CGAL::Triangulation_data_structure_2<Vb, Fb>;
  using Delaunay    = CGAL::Delaunay_triangulation_2<Kernel, Tds>;
  using Alpha_shape = CGAL::Alpha_shape_2<Delaunay>;

  // using Neighbor_query_3 = CGAL::Shape_detection::Point_set::
  //   Sphere_neighbor_query<Kernel, Indices, Point_map_3>;

  using Neighbor_query_3 = CGAL::Shape_detection::Point_set::
    K_neighbor_query<Kernel, Indices, Point_map_3>;
  using Planar_region    = CGAL::Shape_detection::Point_set::
    Least_squares_plane_fit_region<Kernel, Indices, Point_map_3, Vector_map_3>;
  using Planar_sorting   = CGAL::Shape_detection::Point_set::
    Least_squares_plane_fit_sorting<Kernel, Indices, Neighbor_query_3, Point_map_3>;
  using Region_growing_3 = CGAL::Shape_detection::
    Region_growing<Indices, Neighbor_query_3, Planar_region, typename Planar_sorting::Seed_map>;

  using Points_2           = std::vector<Point_2>;
  using Pair_item_2        = std::pair<Point_2, Vector_2>;
  using Pair_range_2       = std::vector<Pair_item_2>;
  using First_of_pair_map  = CGAL::First_of_pair_property_map<Pair_item_2>;
  using Second_of_pair_map = CGAL::Second_of_pair_property_map<Pair_item_2>;

  using Identity_map_2   = CGAL::
    Identity_property_map<Point_2>;
  using Neighbor_query_2 = CGAL::Shape_detection::Point_set::
    Sphere_neighbor_query<Kernel, Points_2, Identity_map_2>;

  // using Neighbor_query_2 = CGAL::Shape_detection::Point_set::
  //   K_neighbor_query<Kernel, Points_2, Identity_map_2>;

  using Estimate_normals_2 = KSR::
    Estimate_normals_2<Kernel, Points_2, Neighbor_query_2>;
  using Linear_region      = CGAL::Shape_detection::Point_set::
    Least_squares_line_fit_region<Kernel, Pair_range_2, First_of_pair_map, Second_of_pair_map>;
  using Linear_sorting     = CGAL::Shape_detection::Point_set::
    Least_squares_line_fit_sorting<Kernel, Points_2, Neighbor_query_2, Identity_map_2>;
  using Region_growing_2   = CGAL::Shape_detection::
    Region_growing<Points_2, Neighbor_query_2, Linear_region, typename Linear_sorting::Seed_map>;

  using Visibility_label = KSR::Visibility_label;
  using Visibility       = KSR_3::Visibility<Kernel, Point_map_3, Vector_map_3>;
  using Graphcut         = KSR_3::Graphcut<Kernel>;

public:

  Reconstruction(
    const Input_range& input_range,
    const Point_map& point_map,
    const Vector_map& normal_map,
    const Semantic_map& semantic_map,
    Data_structure& data,
    const bool verbose,
    const bool debug = true) :
  m_input_range(input_range),
  m_point_map(point_map),
  m_normal_map(normal_map),
  m_semantic_map(semantic_map),
  m_point_map_3(m_input_range, m_point_map),
  m_normal_map_3(m_input_range, m_normal_map),
  m_data(data),
  m_debug(debug),
  m_verbose(verbose),
  m_planar_shape_type(Planar_shape_type::CONVEX_HULL) {

    clear();
    collect_points(Semantic_label::GROUND           , m_ground_points);
    collect_points(Semantic_label::BUILDING_BOUNDARY, m_boundary_points);
    collect_points(Semantic_label::BUILDING_INTERIOR, m_interior_points);

    bool is_ground = (m_ground_points.size() >= 3);
    bool is_wall = (m_boundary_points.size() >= 3);
    bool is_roof = (m_interior_points.size() >= 3);

    if (!is_ground && !is_wall && !is_roof) {
      collect_all_points(m_free_form_points);
      if (m_verbose) {
        std::cout << std::endl << "--- FREE-FORM RECONSTRUCTION: " << std::endl;
        std::cout << "* num points: " << m_free_form_points.size() << std::endl;
      }
      return;
    }

    if (!is_roof) {
      CGAL_assertion_msg(false, "TODO: ADD NEW RECONSTRUCTION TYPE! WALLS / GROUND / TREES?");
    }
    CGAL_assertion(is_roof);

    if (!is_ground) {
      if (is_wall) {
        get_ground_points_from_walls();
      } else {
        get_ground_points_from_roofs();
      }
      is_ground = true;
    }
    CGAL_assertion(is_ground);

    if (m_verbose) {
      std::cout << std::endl << "--- BUILDING RECONSTRUCTION: " << std::endl;
      std::cout << "* num ground points: "   << m_ground_points.size()   << std::endl;
      std::cout << "* num boundary points: " << m_boundary_points.size() << std::endl;
      std::cout << "* num interior points: " << m_interior_points.size() << std::endl;
    }
  }

  template<typename NamedParameters>
  bool detect_planar_shapes(
    const std::string& file_name,
    const NamedParameters& np) {

    if (m_verbose) {
      std::cout << std::endl << "--- DETECTING PLANAR SHAPES: " << std::endl;
    }
    m_planes.clear();
    m_polygons.clear();
    m_region_map.clear();

    if (m_free_form_points.size() == 0) {
      create_ground_plane();
      create_approximate_walls(np);
      create_approximate_roofs(np);
    } else {
      create_all_planar_shapes(np, file_name);
    }

    CGAL_assertion(m_planes.size() == m_polygons.size());
    CGAL_assertion(m_polygons.size() == m_region_map.size());
    //if (m_debug)
    dump_polygons("detected-planar-shapes");

    if (m_polygons.size() == 0) {
      if (m_verbose) std::cout << "* no planar shapes found" << std::endl;
      return false;
    }
    return true;
  }

  template<typename RegionMap, typename NamedParameters>
  bool planar_shapes_from_map(
    RegionMap &region_map,
    const NamedParameters& np
  ) {

    if (m_verbose) {
      std::cout << std::endl << "--- PLANAR SHAPES from map: " << std::endl;
    }
    m_planes.clear();
    m_polygons.clear();
    m_region_map.clear();

    std::vector<std::vector<std::size_t> > regions;

    for (std::size_t i = 0; i < m_input_range.size(); i++) {
     if (region_map[i] < 0)
       continue;
      if (regions.size() <= region_map[i])
        regions.resize(region_map[i] + 1);
      regions[region_map[i]].push_back(i);
    }

    for (const auto& region : regions) {
      const auto plane = fit_plane(region);
      const std::size_t shape_idx = add_planar_shape(region, plane);
      CGAL_assertion(shape_idx != std::size_t(-1));
      m_region_map[shape_idx] = region;
    }

    return true;
  }

  template<typename NamedParameters>
  bool regularize_planar_shapes(
    const NamedParameters& np) {

    if (m_verbose) {
      std::cout << std::endl << "--- REGULARIZING PLANAR SHAPES: " << std::endl;
    }

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
    const FT max_accepted_angle    = FT(10); // parameters::choose_parameter(
      // parameters::get_parameter(np, internal_np::angle_threshold), FT(15));
    const FT max_distance_to_plane = FT(1) / FT(5); // parameters::choose_parameter(
      // parameters::get_parameter(np, internal_np::distance_threshold), FT(1));
    const Vector_3 symmetry_axis(FT(0), FT(0), FT(1));

    std::vector<Plane_3> planes;
    std::vector<Indices> regions;
    create_planes_and_regions(planes, regions);

    CGAL_assertion(planes.size() > 0);
    CGAL_assertion(planes.size() == regions.size());

    Plane_map plane_map;
    Point_to_plane_map point_to_plane_map(m_input_range, regions);
    CGAL::Shape_regularization::Planes::regularize_planes(
      m_input_range,
      m_point_map,
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
      const auto& plane  = planes[i];
      const auto& region = regions[i];

      const std::size_t shape_idx = add_planar_shape(region, plane);
      CGAL_assertion(shape_idx != std::size_t(-1));
      m_region_map[shape_idx] = region;
    }
    CGAL_assertion(m_polygons.size() == num_polygons);
    CGAL_assertion(m_polygons.size() == m_planes.size());
    CGAL_assertion(m_polygons.size() == m_region_map.size());

    if (m_verbose) {
      std::cout << "* num regularized planes: " << m_planes.size() << std::endl;
    }

    if (m_debug) dump_polygons("regularized-planar-shapes");
    return true;
  }

  template<typename NamedParameters>
  bool compute_model(
    const NamedParameters& np) {

    if (m_verbose) {
      std::cout << std::endl << "--- COMPUTING THE MODEL: " << std::endl;
    }

    if (m_data.number_of_volumes() == 0) {
      if (m_verbose) std::cout << "* no volumes found, skipping" << std::endl;
      return false;
    }

    if (m_verbose) std::cout << "* computing visibility ... ";
    std::map<PFace, Indices> pface_points;
    assign_points_to_pfaces(pface_points);
    const Visibility visibility(
      m_data, pface_points, m_point_map_3, m_normal_map_3);

    CGAL_assertion(m_data.volumes().size() > 0);
    visibility.compute(m_data.volumes());
    //dump_visibility("visibility/visibility", pface_points);

    if (m_verbose) {
      std::cout << "done" << std::endl;
      std::cout << "* applying graphcut ... ";
    }

    const FT beta = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::graphcut_beta), FT(1) / FT(2));

    Graphcut graphcut(m_data, beta);
    graphcut.compute(m_data.volumes(), visibility.inliers());
    //dump_volumes("graphcut/graphcut");

    if (m_verbose) {
      std::cout << "done" << std::endl;
      std::cout << "* extracting the model ... ";
    }

    extract_surface_model();
    if (m_debug) dump_model("reconstructed-model");

    if (m_verbose) std::cout << "done" << std::endl;
    return true;
  }

  const std::vector<Polygon_3>& planar_shapes() const {
    return m_polygons;
  }

  const Polygon_map& polygon_map() const {
    return m_polygon_map;
  }

  void clear() {
    m_ground_points.clear();
    m_boundary_points.clear();
    m_interior_points.clear();
    m_free_form_points.clear();
    m_region_map.clear();
    m_polygons.clear();
    m_planes.clear();
  }

private:
  const Input_range& m_input_range;
  const Point_map& m_point_map;
  const Vector_map& m_normal_map;
  const Semantic_map& m_semantic_map;

  Point_map_3  m_point_map_3;
  Vector_map_3 m_normal_map_3;

  Data_structure& m_data;
  const bool m_debug;
  const bool m_verbose;
  const Planar_shape_type m_planar_shape_type;
  const Converter m_converter;

  std::vector<std::size_t> m_ground_points;
  std::vector<std::size_t> m_boundary_points;
  std::vector<std::size_t> m_interior_points;
  std::vector<std::size_t> m_free_form_points;

  std::vector<Polygon_3> m_polygons;
  std::vector<Plane_3> m_planes;
  Polygon_map m_polygon_map;

  std::map<std::size_t, Indices> m_region_map;

  void collect_points(
    const Semantic_label output_label,
    std::vector<std::size_t>& indices) const {

    indices.clear();
    for (std::size_t i = 0; i < m_input_range.size(); ++i) {
      const Semantic_label label =
      get(m_semantic_map, *(m_input_range.begin() + i));
      if (label == output_label)
        indices.push_back(i);
    }
  }

  void collect_all_points(std::vector<std::size_t>& indices) const {

    indices.clear();
    indices.reserve(m_input_range.size());
    for (std::size_t i = 0; i < m_input_range.size(); ++i) {
      indices.push_back(i);
    }
  }

  void get_ground_points_from_walls() {

    CGAL_assertion(m_ground_points.size()    < 3);
    CGAL_assertion(m_boundary_points.size() >= 3);
    if (m_verbose) std::cout << "* getting ground points from facade points" << std::endl;
    get_zero_level_points(m_boundary_points, m_ground_points);
    CGAL_assertion(m_ground_points.size() >= 3);

    // CGAL_assertion_msg(false, "TODO: ADD MISSING GROUND POINTS, GET FROM WALLS!");
  }

  void get_ground_points_from_roofs() {

    CGAL_assertion(m_ground_points.size()    < 3);
    CGAL_assertion(m_interior_points.size() >= 3);
    if (m_verbose) std::cout << "* getting ground points from roof points" << std::endl;
    get_zero_level_points(m_interior_points, m_ground_points);
    CGAL_assertion(m_ground_points.size() >= 3);

    // CGAL_assertion_msg(false, "TODO: ADD MISSING GROUND POINTS, GET FROM ROOFS!");
  }

  void get_zero_level_points(
    const std::vector<std::size_t>& input_range,
    std::vector<std::size_t>& output) const {

    CGAL_assertion(input_range.size() >= 3);
    output.clear();

    FT min_z = +FT(1000000000000);
    FT max_z = -FT(1000000000000);
    for (const std::size_t idx : input_range) {
      CGAL_assertion(idx < m_input_range.size());
      const auto& point = get(m_point_map_3, idx);
      min_z = CGAL::min(min_z, point.z());
      max_z = CGAL::max(max_z, point.z());
    }
    CGAL_assertion(min_z < +FT(1000000000000));
    CGAL_assertion(max_z > -FT(1000000000000));
    CGAL_assertion(max_z > min_z);

    const FT d = (max_z - min_z) / FT(100);
    const FT top_level = min_z + d;

    for (const std::size_t idx : input_range) {
      CGAL_assertion(idx < m_input_range.size());
      const auto& point = get(m_point_map_3, idx);
      if (point.z() < top_level) output.push_back(idx);
    }
    CGAL_assertion(output.size() >= 3);
  }

  void create_ground_plane() {

    if (m_verbose) std::cout << "* creating ground plane ... ";
    if (m_ground_points.size() < 3) {
      if (m_verbose) std::cout << "omitted, no points available" << std::endl;
      return;
    }

    const auto plane = fit_plane(m_ground_points);
    const std::size_t shape_idx = add_planar_shape(m_ground_points, plane);
    CGAL_assertion(shape_idx != std::size_t(-1));
    m_region_map[shape_idx] = m_ground_points;
    extend_ground_plane(shape_idx);
    if (m_verbose) std::cout << "done" << std::endl;
  }

  void extend_ground_plane(const std::size_t shape_idx) {

    FT min_x = +FT(1000000000000), min_y = +FT(1000000000000);
    FT max_x = -FT(1000000000000), max_y = -FT(1000000000000);
    CGAL_assertion(m_interior_points.size() >= 3);
    for (const std::size_t idx : m_interior_points) {
      CGAL_assertion(idx < m_input_range.size());
      const auto& point = get(m_point_map_3, idx);
      min_x = CGAL::min(min_x, point.x());
      min_y = CGAL::min(min_y, point.y());
      max_x = CGAL::max(max_x, point.x());
      max_y = CGAL::max(max_y, point.y());
    }

    const Point_3 a(min_x, min_y, FT(0));
    const Point_3 b(max_x, min_y, FT(0));
    const Point_3 c(max_x, max_y, FT(0));
    const Point_3 d(min_x, max_y, FT(0));

    const auto& plane = m_planes[shape_idx];
    const auto p0 = plane.projection(a);
    const auto p1 = plane.projection(b);
    const auto p2 = plane.projection(c);
    const auto p3 = plane.projection(d);

    auto& polygon = m_polygons[shape_idx];
    polygon.clear();
    polygon.push_back(p0);
    polygon.push_back(p1);
    polygon.push_back(p2);
    polygon.push_back(p3);
  }

  const Plane_3 fit_plane(const std::vector<std::size_t>& region) const {

    std::vector<IPoint_3> points;
    points.reserve(region.size());
    for (const std::size_t idx : region) {
      CGAL_assertion(idx < m_input_range.size());
      points.push_back(m_converter(get(m_point_map_3, idx)));
    }
    CGAL_assertion(points.size() == region.size());

    IPlane_3 fitted_plane;
    IPoint_3 fitted_centroid;
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

  std::size_t add_planar_shape(
    const std::vector<std::size_t>& region, const Plane_3& plane) {

    switch (m_planar_shape_type) {
      case Planar_shape_type::CONVEX_HULL: {
        return add_convex_hull_shape(region, plane);
      }
      case Planar_shape_type::RECTANGLE: {
        return add_rectangle_shape(region, plane);
      }
      default: {
        CGAL_assertion_msg(false, "ERROR: ADD PLANAR SHAPE, WRONG TYPE!");
        return std::size_t(-1);
      }
    }
    return std::size_t(-1);
  }

  std::size_t add_convex_hull_shape(
    const std::vector<std::size_t>& region, const Plane_3& plane) {

    std::vector<Point_2> points;
    points.reserve(region.size());
    for (const std::size_t idx : region) {
      CGAL_assertion(idx < m_input_range.size());
      const auto& p = get(m_point_map_3, idx);
      const auto q = plane.projection(p);
      const auto point = plane.to_2d(q);
      points.push_back(point);
    }
    CGAL_assertion(points.size() == region.size());

    std::vector<Point_2> ch;
    CGAL::convex_hull_2(points.begin(), points.end(), std::back_inserter(ch) );

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

  std::size_t add_rectangle_shape(
    const std::vector<std::size_t>& /* region */, const Plane_3& /* plane */) {

    CGAL_assertion_msg(false, "TODO: ADD RECTANGLE SHAPE!");
    return std::size_t(-1);
  }

  template<typename NamedParameters>
  void create_all_planar_shapes(const NamedParameters& np, const std::string& file_name) {

    if (m_free_form_points.size() < 3) {
      if (m_verbose) std::cout << "* no points found, skipping" << std::endl;
      return;
    }
    if (m_verbose) std::cout << "* getting planar shapes using region growing" << std::endl;
    const std::size_t num_shapes = compute_planar_shapes_with_rg(np, m_free_form_points, file_name);
    if (m_verbose) std::cout << "* found " << num_shapes << " approximate walls" << std::endl;
  }

  template<typename NamedParameters>
  void create_approximate_walls(const NamedParameters& np) {

    if (m_boundary_points.size() < 3) {
      create_walls_from_roof_boundaries(np);
      return;
    }
    if (m_verbose) std::cout << "* getting walls using region growing" << std::endl;
    const std::size_t num_shapes = compute_planar_shapes_with_rg(np, m_boundary_points, "");
    if (m_verbose) std::cout << "* found " << num_shapes << " approximate walls" << std::endl;
  }

  template<typename NamedParameters>
  void create_approximate_roofs(const NamedParameters& np) {

    if (m_interior_points.size() < 3) {
      if (m_verbose) std::cout << "* no roof points found, skipping" << std::endl;
      return;
    }
    if (m_verbose) std::cout << "* getting roofs using region growing" << std::endl;
    const std::size_t num_shapes = compute_planar_shapes_with_rg(np, m_interior_points, "");
    if (m_verbose) std::cout << "* found " << num_shapes << " approximate roofs" << std::endl;
  }

  template<typename NamedParameters>
  std::size_t compute_planar_shapes_with_rg(
    const NamedParameters& np,
    const std::vector<std::size_t>& input_range,
    const std::string& file_name) {

    std::vector< std::vector<std::size_t> > regions;
    apply_region_growing_3(np, input_range, regions);

    if (file_name.size() > 0)
      dump(m_input_range, m_point_map_3, m_normal_map_3, regions, file_name);

    for (const auto& region : regions) {
      const auto plane = fit_plane(region);
      const std::size_t shape_idx = add_planar_shape(region, plane);
      CGAL_assertion(shape_idx != std::size_t(-1));
      m_region_map[shape_idx] = region;
    }
    return regions.size();
  }

  template<typename NamedParameters>
  void apply_region_growing_3(
    const NamedParameters& np,
    const std::vector<std::size_t>& input_range,
    std::vector< std::vector<std::size_t> >& regions) const {

    // const FT radius = parameters::choose_parameter(
    //   parameters::get_parameter(np, internal_np::neighbor_radius), FT(1));

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
    Neighbor_query_3 neighbor_query(input_range, k, m_point_map_3);

    Planar_region planar_region(input_range,
      max_distance_to_plane, max_accepted_angle, min_region_size,
      m_point_map_3, m_normal_map_3);

    Planar_sorting sorting(
      input_range, neighbor_query, m_point_map_3);
    sorting.sort();

    std::vector<Indices> result;
    Region_growing_3 region_growing(
      input_range, neighbor_query, planar_region, sorting.seed_map());
    region_growing.detect(std::back_inserter(result));

    // Convert indices.
    regions.clear();
    regions.reserve(result.size());

    Indices region;
    for (const auto& indices : result) {
      region.clear();
      for (const std::size_t index : indices) {
        region.push_back(input_range[index]);
      }
      regions.push_back(region);
    }
    CGAL_assertion(regions.size() == result.size());
  }

  template<typename NamedParameters>
  void create_walls_from_roof_boundaries(const NamedParameters& np) {

    if (m_interior_points.size() < 3) {
      if (m_verbose) std::cout << "* no facade points found, skipping" << std::endl;
      return;
    }

    CGAL_assertion(m_interior_points.size() >= 3);
    if (m_verbose) std::cout << "* getting walls using roof boundaries" << std::endl;

    const FT max_accepted_angle = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::angle_threshold), FT(15));
    std::vector<std::size_t> wall_points, roof_points;
    split_points(max_accepted_angle, m_interior_points, wall_points, roof_points);
    // dump_points(wall_points, "wall-points");
    // dump_points(roof_points, "roof-points");

    std::size_t num_shapes = 0;
    if (wall_points.size() >= 3) {
      num_shapes += compute_planar_shapes_with_rg(np, wall_points, "");
      // dump_polygons("walls-1");
    }

    if (roof_points.size() >= 3) {
      num_shapes += add_polygons_using_alpha_shapes(np, roof_points);
      // dump_polygons("walls-2");
    }

    if (m_verbose) std::cout << "* found " << num_shapes << " approximate walls" << std::endl;
    // CGAL_assertion_msg(false, "TODO: GET WALLS FROM ROOF BOUNDARIES!");
  }

  void split_points(
    const FT max_accepted_angle,
    const std::vector<std::size_t>& all_points,
    std::vector<std::size_t>& wall_points,
    std::vector<std::size_t>& roof_points) const {

    wall_points.clear(); roof_points.clear();
    const Vector_3 ref = Vector_3(FT(0), FT(0), FT(1));
    for (const std::size_t idx : all_points) {
      CGAL_assertion(idx < m_input_range.size());
      const auto& normal = get(m_normal_map_3, idx);

      FT angle = approximate_angle(normal, ref);
      if (angle > FT(90)) angle = FT(180) - angle;
      angle = FT(90) - angle;
      if (angle <= max_accepted_angle) wall_points.push_back(idx);
      else roof_points.push_back(idx);
    }
  }

  template<typename NamedParameters>
  std::size_t add_polygons_using_alpha_shapes(
    const NamedParameters& np,
    const std::vector<std::size_t>& input_range) {

    Delaunay triangulation;
    CGAL_assertion(input_range.size() >= 3);
    create_triangulation(input_range, triangulation);
    if (triangulation.number_of_faces() == 0) return 0;

    std::vector<Point_2> boundary_points;
    add_filtered_points(np, triangulation, boundary_points);
    // dump_points(boundary_points, "boundary-points");

    std::vector<Indices> regions;
    apply_region_growing_2(np, boundary_points, regions);
    // dump_points(boundary_points, regions, "boundary-regions");

    std::vector<Line_2> lines;
    create_lines(boundary_points, regions, lines);
    CGAL_assertion(lines.size() == regions.size());
    // dump_points(boundary_points, regions, lines, "projected-regions");

    std::vector<Segment_2> segments;
    create_segments(boundary_points, regions, lines, segments);
    // dump_segments(segments, "boundary-segments");

    const std::size_t num_walls = add_walls_from_segments(segments);
    // std::cout << "num walls: " << num_walls << std::endl;
    return num_walls;
  }

  void create_triangulation(
    const std::vector<std::size_t>& input_range,
    Delaunay& triangulation) {

    triangulation.clear();
    for (const std::size_t idx : input_range) {
      CGAL_assertion(idx < m_input_range.size());
      const auto& point = get(m_point_map_3, idx);
      const auto vh = triangulation.insert(
        KSR::point_2_from_point_3(point));
      vh->info().z = point.z();
    }
  }

  template<typename NamedParameters>
  void add_filtered_points(
    const NamedParameters& np,
    Delaunay& triangulation,
    std::vector<Point_2>& boundary_points) {

    CGAL_precondition(triangulation.number_of_faces() != 0);
    const FT distance_threshold = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::distance_threshold), FT(1));
    const FT alpha = distance_threshold / FT(2);
    CGAL_precondition(alpha > FT(0));

    Alpha_shape alpha_shape(triangulation, alpha, Alpha_shape::GENERAL);
    sample_edges(np, alpha_shape, boundary_points);
  }

  template<typename NamedParameters>
  void sample_edges(
    const NamedParameters& np,
    const Alpha_shape& alpha_shape,
    std::vector<Point_2>& boundary_points) const {

    const FT distance_threshold = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::distance_threshold), FT(1));
    const FT edge_sampling = distance_threshold / FT(4);
    CGAL_precondition(edge_sampling > FT(0));

    for (auto eit = alpha_shape.alpha_shape_edges_begin();
      eit != alpha_shape.alpha_shape_edges_end(); ++eit) {

      const auto& source = eit->first->vertex((eit->second + 1) % 3)->point();
      const auto& target = eit->first->vertex((eit->second + 2) % 3)->point();
      sample_edge(edge_sampling, source, target, boundary_points);
    }
  }

  void sample_edge(
    const FT edge_sampling,
    const Point_2& source, const Point_2& target,
    std::vector<Point_2>& boundary_points) const {

    CGAL_precondition(edge_sampling > FT(0));
    const FT distance = KSR::distance(source, target);
    const std::size_t nb_pts = static_cast<std::size_t>(
      CGAL::to_double(distance / edge_sampling)) + 1;

    CGAL_precondition(nb_pts > 0);
    for (std::size_t i = 0; i <= nb_pts; ++i) {
      const FT ratio = static_cast<FT>(i) / static_cast<FT>(nb_pts);
      boundary_points.push_back(
        Point_2(
          source.x() * (FT(1) - ratio) + target.x() * ratio,
          source.y() * (FT(1) - ratio) + target.y() * ratio));
    }
  }

  template<typename NamedParameters>
  void apply_region_growing_2(
    const NamedParameters& np,
    const std::vector<Point_2>& input_range,
    std::vector< std::vector<std::size_t> >& regions) const {

    const FT distance_threshold = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::distance_threshold), FT(1));
    CGAL_precondition(distance_threshold > FT(0));
    const FT angle_threshold = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::angle_threshold), FT(15));
    CGAL_precondition(angle_threshold > FT(0));
    const std::size_t min_region_size = 20;
    CGAL_precondition(min_region_size > 0);

    regions.clear();
    Identity_map_2 identity_map_2;
    const FT scale = distance_threshold * FT(2);
    Neighbor_query_2 neighbor_query(input_range, scale, identity_map_2);

    std::vector<Vector_2> normals;
    Estimate_normals_2 estimator(input_range, neighbor_query);
    estimator.get_normals(normals);
    CGAL_assertion(input_range.size() == normals.size());

    Pair_range_2 range;
    range.reserve(input_range.size());
    for (std::size_t i = 0; i < input_range.size(); ++i)
      range.push_back(std::make_pair(input_range[i], normals[i]));

    First_of_pair_map point_map;
    Second_of_pair_map normal_map;
    Linear_region region(
      range, distance_threshold, angle_threshold, min_region_size,
      point_map, normal_map);

    Linear_sorting sorting(
      input_range, neighbor_query, identity_map_2);
    sorting.sort();

    Region_growing_2 region_growing(
      input_range, neighbor_query, region, sorting.seed_map());
    region_growing.detect(std::back_inserter(regions));
  }

  void create_lines(
    const std::vector<Point_2>& input_range,
    const std::vector< std::vector<std::size_t> >& regions,
    std::vector<Line_2>& lines) const {

    lines.clear();
    lines.reserve(regions.size());
    for (const auto& region : regions) {
      const auto line = fit_line(input_range, region);
      lines.push_back(line);
    }
    CGAL_assertion(lines.size() == regions.size());
  }

  const Line_2 fit_line(
    const std::vector<Point_2>& input_range,
    const std::vector<std::size_t>& region) const {

    std::vector<IPoint_2> points;
    points.reserve(region.size());
    for (const std::size_t idx : region) {
      CGAL_assertion(idx < input_range.size());
      points.push_back(m_converter(input_range[idx]));
    }
    CGAL_assertion(points.size() == region.size());

    ILine_2 fitted_line;
    IPoint_2 fitted_centroid;
    CGAL::linear_least_squares_fitting_2(
      points.begin(), points.end(),
      fitted_line, fitted_centroid,
      CGAL::Dimension_tag<0>());

    const Line_2 line(
      static_cast<FT>(fitted_line.a()),
      static_cast<FT>(fitted_line.b()),
      static_cast<FT>(fitted_line.c()));
    return line;
  }

  void create_segments(
    const std::vector<Point_2>& input_range,
    const std::vector< std::vector<std::size_t> >& regions,
    const std::vector<Line_2>& lines,
    std::vector<Segment_2>& segments) const {

    CGAL_assertion(lines.size() == regions.size());
    CGAL_assertion(m_planes.size() > 0);

    segments.clear();
    segments.reserve(lines.size());
    for (std::size_t i = 0; i < lines.size(); ++i) {
      Point_2 source, target;
      KSR::boundary_points_on_line_2(
        input_range, regions[i], lines[i], source, target);
      segments.push_back(Segment_2(source, target));
    }
    CGAL_assertion(segments.size() == lines.size());
  }

  std::size_t add_walls_from_segments(
    const std::vector<Segment_2>& segments) {

    FT min_z = +FT(1000000000000);
    FT max_z = -FT(1000000000000);

    for (const std::size_t idx : m_boundary_points) {
      CGAL_assertion(idx < m_input_range.size());
      const auto& point = get(m_point_map_3, idx);
      min_z = CGAL::min(min_z, point.z());
      max_z = CGAL::max(max_z, point.z());
    }

    for (const std::size_t idx : m_interior_points) {
      CGAL_assertion(idx < m_input_range.size());
      const auto& point = get(m_point_map_3, idx);
      min_z = CGAL::min(min_z, point.z());
      max_z = CGAL::max(max_z, point.z());
    }

    CGAL_assertion(min_z <= max_z);
    for (const auto& segment : segments) {
      const auto& source = segment.source();
      const auto& target = segment.target();

      const Point_3 a(source.x(), source.y(), min_z);
      const Point_3 b(target.x(), target.y(), min_z);
      const Point_3 c(target.x(), target.y(), max_z);
      const Point_3 d(source.x(), source.y(), max_z);

      const std::size_t shape_idx = m_polygons.size();
      m_polygons.push_back({a, b, c, d});
      m_planes.push_back(Plane_3(a, b, c));
      m_region_map[shape_idx] = std::vector<std::size_t>();
    }

    return segments.size();
  }

  void create_planes_and_regions(
    std::vector<Plane_3>& planes,
    std::vector<Indices>& regions) const {

    planes.clear();
    planes.reserve(m_region_map.size());

    regions.clear();
    regions.reserve(m_region_map.size());

    CGAL_assertion(m_planes.size() == m_region_map.size());
    for (const auto& item : m_region_map) {
      const std::size_t shape_idx = item.first;

      const auto& plane = m_planes[shape_idx];
      CGAL_assertion(plane != Plane_3());
      planes.push_back(plane);

      const auto& region = item.second;
      CGAL_assertion(region.size() > 0);
      regions.push_back(region);
    }
    CGAL_assertion(planes.size()  == m_region_map.size());
    CGAL_assertion(regions.size() == m_region_map.size());
  }

  void assign_points_to_pfaces(std::map<PFace, Indices>& pface_points) const {

    pface_points.clear();
    for (std::size_t i = 0; i < m_data.number_of_support_planes(); ++i) {
      const auto pfaces = m_data.pfaces(i);
      for (const auto pface : pfaces) {
        pface_points[pface] = Indices();
      }
    }

    CGAL_assertion(m_region_map.size() > 0);
    for (const auto& item : m_region_map) {
      const std::size_t shape_idx = item.first;
      const auto& indices = item.second;

      const int sp_idx = m_data.support_plane_index(shape_idx);
      CGAL_assertion(sp_idx >= 6);
      const std::size_t support_plane_idx = static_cast<std::size_t>(sp_idx);
      // dump_points(indices, "sp-points-" + std::to_string(support_plane_idx));

      const auto pfaces = m_data.pfaces(support_plane_idx);
      for (const auto pface : pfaces) {
        const auto pvertices = m_data.pvertices_of_pface(pface);

        Delaunay tri;
        for (const auto pvertex : pvertices) {
          CGAL_assertion(m_data.has_ivertex(pvertex));
          const auto ivertex = m_data.ivertex(pvertex);
          const auto& point = m_data.point_2(support_plane_idx, ivertex);
          tri.insert(point);
        }

        for (const std::size_t index : indices) {
          const auto& point = get(m_point_map_3, index);
          const auto query = m_data.to_2d(support_plane_idx, point);
          const auto fh = tri.locate(query);
          if (fh != nullptr && !tri.is_infinite(fh)) {
            pface_points[pface].push_back(index);
          }
        }
      }
    }

    // for (const auto& item : pface_points) {
    //   dump_points(item.second, "pf-points-" + m_data.str(item.first));
    // }
  }

  void extract_surface_model() {
    create_surface_model();
    orient_surface_model();
  }

  void create_surface_model() {
    auto& model = m_data.reconstructed_model();
    model.clear();

    const auto& volumes = m_data.volumes();
    const auto& items   = m_data.pface_neighbors();

    for (const auto& item : items) {
      const auto& pface = item.first;
      const auto& neighbors = item.second;

      const int idx1 = neighbors.first;
      const int idx2 = neighbors.second;

      // std::cout << "idx1/2: " << idx1 << "/" << idx2 << std::endl;

      CGAL_assertion(idx1 >= 0 || idx2 >= 0);
      if (idx1 >= 0 && idx2 >= 0) {
        const auto& volume1 = volumes[idx1];
        const auto& volume2 = volumes[idx2];

        const auto label1 = volume1.visibility;
        const auto label2 = volume2.visibility;

        if (
          label1 == Visibility_label::INSIDE &&
          label2 == Visibility_label::OUTSIDE) {
          model.pfaces.push_back(pface);
        } else if (
          label1 == Visibility_label::OUTSIDE &&
          label2 == Visibility_label::INSIDE) {
          model.pfaces.push_back(pface);
        }
        continue;
      }

      if (idx1 >= 0) {
        CGAL_assertion(idx2 < 0);
        const auto& volume1 = volumes[idx1];
        const auto label1 = volume1.visibility;
        if (label1 == Visibility_label::INSIDE) {
          model.pfaces.push_back(pface);
        }
        continue;
      }

      if (idx2 >= 0) {
        CGAL_assertion(idx1 < 0);
        const auto& volume2 = volumes[idx2];
        const auto label2 = volume2.visibility;
        if (label2 == Visibility_label::INSIDE) {
          model.pfaces.push_back(pface);
        }
        continue;
      }
    }
  }

  void orient_surface_model() {
    return;
    CGAL_assertion_msg(false, "TODO: ORIENT SURFACE MODEL!");
  }

  void dump_points(
    const std::vector<Point_2>& points,
    const std::string file_name) const {

    KSR_3::Saver<Kernel> saver;
    saver.export_points_2(points, file_name);
  }

  void dump_points(
    const std::vector<Point_2>& boundary_points,
    const std::vector< std::vector<std::size_t> >& regions,
    const std::string file_name) const {

    std::vector<Point_2> points;
    std::vector< std::vector<Point_2> > all_points;
    all_points.reserve(regions.size());

    for (const auto& region : regions) {
      points.clear();
      for (const std::size_t index : region) {
        CGAL_assertion(index < boundary_points.size());
        const auto& point = boundary_points[index];
        points.push_back(point);
      }
      CGAL_assertion(points.size() == region.size());
      all_points.push_back(points);
    }
    CGAL_assertion(all_points.size() == regions.size());

    KSR_3::Saver<Kernel> saver;
    saver.export_points_2(all_points, file_name);
  }

  void dump_points(
    const std::vector<Point_2>& boundary_points,
    const std::vector< std::vector<std::size_t> >& regions,
    const std::vector<Line_2>& lines,
    const std::string file_name) const {

    std::vector<Point_2> points;
    std::vector< std::vector<Point_2> > all_points;
    all_points.reserve(regions.size());

    for (std::size_t i = 0; i < regions.size(); ++i) {
      points.clear();
      for (const std::size_t index : regions[i]) {
        CGAL_assertion(index < boundary_points.size());
        const auto& point = boundary_points[index];
        const auto proj = lines[i].projection(point);
        points.push_back(proj);
      }
      CGAL_assertion(points.size() == regions[i].size());
      all_points.push_back(points);
    }
    CGAL_assertion(all_points.size() == regions.size());

    KSR_3::Saver<Kernel> saver;
    saver.export_points_2(all_points, file_name);
  }

  void dump_points(
    const std::vector<std::size_t>& indices,
    const std::string file_name) const {

    std::vector<Point_3> points;
    points.reserve(indices.size());
    for (const std::size_t index : indices) {
      const auto& point = get(m_point_map_3, index);
      points.push_back(point);
    }
    CGAL_assertion(points.size() == indices.size());

    KSR_3::Saver<Kernel> saver;
    saver.export_points_3(points, file_name);
  }

  void dump_segments(
    const std::vector<Segment_2>& segments,
    const std::string file_name) {

    KSR_3::Saver<Kernel> saver;
    saver.export_segments_2(segments, file_name);
  }

  void dump_polygons(const std::string file_name) {

    KSR_3::Saver<Kernel> saver;
    saver.export_polygon_soup_3(m_polygons, file_name);
  }

  void dump_volumes(const std::string file_name) {

    for (const auto& volume : m_data.volumes()) {
      if (volume.visibility == Visibility_label::INSIDE) {
        dump_volume(m_data, volume.pfaces,
        file_name + "-" + std::to_string(volume.index), false);
      }
    }
  }

  void dump_visibility(const std::string file_name, const std::map<PFace, Indices> &pface_points) {
    for (const auto& volume : m_data.volumes()) {
      std::size_t sample_count = 0;
      for (auto pface : volume.pfaces) {
        const auto indices = pface_points.at(pface);
        sample_count += indices.size();
      }
      dump_visi<Data_structure, PFace, FT>(m_data, volume.pfaces,
        file_name + "-" + std::to_string(volume.index) + "-" + std::to_string(volume.inside_count) + "-"
        + std::to_string(volume.outside_count) + "-" + std::to_string(sample_count), (volume.inside_count)/(volume.inside_count + volume.outside_count));
    }
  }

  void dump_model(const std::string file_name) {

    From_EK from_EK;

    std::vector<Point_3> polygon;
    std::vector< std::vector<Point_3> > polygons;
    const auto& model = m_data.reconstructed_model();
    std::vector<CGAL::Color> colors;

    std::size_t polygon_id = 0;
    KSR_3::Saver<Kernel> saver;
    for (const auto& pface : model.pfaces) {
      const auto pvertices = m_data.pvertices_of_pface(pface);
      polygon.clear();
      for (const auto pvertex : pvertices) {
        CGAL_assertion(m_data.has_ivertex(pvertex));
        const auto ivertex = m_data.ivertex(pvertex);
        const auto point = from_EK(m_data.point_3(ivertex));
        polygon.push_back(point);
      }
      polygons.push_back(polygon);
      colors.push_back(saver.get_idx_color(pface.first));
      ++polygon_id;
    }
    saver.export_polygon_soup_3(polygons, colors, file_name);
  }
};

#endif //DOXYGEN_RUNNING

} // namespace KSR_3
} // namespace CGAL

#endif // CGAL_KSR_3_RECONSTRUCTION_H
