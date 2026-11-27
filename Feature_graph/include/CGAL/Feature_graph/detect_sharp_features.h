// Copyright (c) 2026 GeometryFactory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Ange Clement

#ifndef CGAL_FEATURE_GRAPH_DETECT_SHARP_FEATURES_H
#define CGAL_FEATURE_GRAPH_DETECT_SHARP_FEATURES_H

#include <CGAL/license/Feature_graph.h>

#include <algorithm>
#include <array>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <type_traits>
#include <variant>
#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/property_map/vector_property_map.hpp>
#include <boost/range/iterator_range.hpp>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <CGAL/assertions.h>
#include <CGAL/number_utils.h>

#include <CGAL/Feature_graph/internal/image_Khalimsky_coordinates.h>
#include <CGAL/Feature_graph/internal/graph_operation.h>
#include <CGAL/Feature_graph/internal/feature_graph_operation.h>
#include <CGAL/Feature_graph/internal/polyline_regularizer.h>
#ifdef CGAL_DGTAL_ENABLED
#include <CGAL/Feature_graph/Image_AmbrosioTortorelli.h>
#endif
#include <CGAL/Feature_graph/internal/convert_image_to_surface_mesh.h>

#include <CGAL/Feature_graph/internal/debugging_methods.h>

namespace CGAL
{
namespace Feature_graph
{
namespace internal
{

template <typename K, typename PolygonMesh, typename Image>
struct Sharpness_estimator_on_image_translator
{
  typedef typename K::Point_3 Point_3;
  typedef typename K::Vector_3 Vector_3;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;

  const Image& image;
  const Point_3 image_origin;
  const Point_3 voxel_size;

  Sharpness_estimator_on_image_translator(const Image& image)
    : image(image), image_origin(image.tx(), image.ty(), image.tz()),
      voxel_size(image.vx(), image.vy(), image.vz())
  { }

  std::size_t get_Khalimsky_index(const Point_3& point) const
  {
    return get_Khalimsky_index_from_image_coordinates(
      image_origin, voxel_size, image.xdim(), image.ydim(), point);
  }

  template <typename PointPropertyMap>
  std::size_t translate_vertex_element(
    const vertex_descriptor& vertex,
    const PolygonMesh& /*pmesh*/,
    const PointPropertyMap& points) const
  {
    return get_Khalimsky_index(get(points, vertex));
  }

  template <typename PointPropertyMap>
  std::size_t translate_halfedge_element(
    const halfedge_descriptor& edge,
    const PolygonMesh& pmesh,
    const PointPropertyMap& points) const
  {
    Point_3 midpoint =
      CGAL::midpoint(get(points, source(edge, pmesh)), get(points, target(edge, pmesh)));
    return get_Khalimsky_index(midpoint);
  }

  template <typename PointPropertyMap>
  std::size_t translate_face_element(
    const face_descriptor& face, const PolygonMesh& pmesh, const PointPropertyMap& points) const
  {
    Vector_3 barycenter(0,0,0);
    int number_of_edges = 0;
    halfedge_descriptor start_he = halfedge(face, pmesh);
    halfedge_descriptor he = start_he;
    do
    {
      barycenter += Vector_3(CGAL::ORIGIN, get(points, source(he, pmesh)));
      number_of_edges++;
      he = next(he, pmesh);
    } while (he != start_he);
    if (number_of_edges != 0)
      barycenter /= number_of_edges;
    return get_Khalimsky_index(CGAL::ORIGIN + barycenter);
  }

  Image translate_domain(const PolygonMesh& /*pmesh*/) const
  {
    return image;
  }
};

template <bool IsthmusThresholdSpecified, bool SimpleThresholdSpecified,
          bool CornerThresholdSpecified, typename BaselineThreshold,
          typename IsthmusThreshold, typename SimpleThreshold, typename CornerThreshold>
struct Line_regularization_parameters
{
  static constexpr bool isthmus_threshold_specified = IsthmusThresholdSpecified;
  static constexpr bool simple_threshold_specified = SimpleThresholdSpecified;
  static constexpr bool corner_threshold_specified = CornerThresholdSpecified;

  BaselineThreshold line_distance_threshold;
  IsthmusThreshold isthmus_line_distance_threshold;
  SimpleThreshold simple_line_distance_threshold;
  CornerThreshold corner_line_distance_threshold;
};

template <typename K, typename PolygonMesh, typename MeshPointPropertyMap,
          typename ElementAndDomainTranslator, typename VertexPointMap,
          typename VertexDimensionMap, typename VertexToMeshElementMap, typename Domain,
          typename RegularizationParameters>
class Line_distance_threshold_provider
{
  using FT = typename K::FT;
  using Vertex_descriptor =
    typename boost::graph_traits<PolygonMesh>::vertex_descriptor;
  using Halfedge_descriptor =
    typename boost::graph_traits<PolygonMesh>::halfedge_descriptor;
  using Face_descriptor =
    typename boost::graph_traits<PolygonMesh>::face_descriptor;
  using Descriptor_variant =
    std::variant<Vertex_descriptor, Halfedge_descriptor, Face_descriptor>;

  const PolygonMesh& pmesh_;
  const MeshPointPropertyMap& mesh_points_;
  const ElementAndDomainTranslator& translator_;
  const VertexPointMap& vertex_points_;
  const VertexDimensionMap& vertex_dimensions_;
  const VertexToMeshElementMap& vertex_to_mesh_element_;
  const Domain& domain_;
  const RegularizationParameters& regularization_parameters_;

  template <typename Field, typename Vertex>
  FT evaluate_field(const Field& field, const Vertex& vertex) const
  {
    const auto& element_holder = get(vertex_to_mesh_element_, vertex);
    const Descriptor_variant& descriptor =
      static_cast<const Descriptor_variant&>(element_holder);
    const char dimension = get(vertex_dimensions_, vertex);
    const auto& point = get(vertex_points_, vertex);
    if (dimension == 0)
      return field.template operator()<CGAL::Dimension_tag<0>>(
        point,
        translator_.translate_vertex_element(std::get<0>(descriptor), pmesh_, mesh_points_),
        domain_);
    if (dimension == 1)
      return field.template operator()<CGAL::Dimension_tag<1>>(
        point,
        translator_.translate_halfedge_element(std::get<1>(descriptor), pmesh_, mesh_points_),
        domain_);
    CGAL_assertion(dimension == 2);
    return field.template operator()<CGAL::Dimension_tag<2>>(
      point,
      translator_.translate_face_element(std::get<2>(descriptor), pmesh_, mesh_points_),
      domain_);
  }

  template <typename Threshold, typename Vertex>
  double evaluate(const Threshold& threshold, const Vertex& vertex) const
  {
    if constexpr (std::is_convertible<Threshold, FT>::value)
      return CGAL::to_double(static_cast<FT>(threshold));
    else
      return CGAL::to_double(evaluate_field(threshold, vertex));
  }

  template <bool ThresholdSpecified, typename Threshold, typename Vertex>
  double effective_threshold(const Threshold& threshold, const Vertex& vertex) const
  {
    double value;
    if constexpr (ThresholdSpecified)
    {
      const double override_value = evaluate(threshold, vertex);
      value = override_value < 0
        ? evaluate(regularization_parameters_.line_distance_threshold, vertex)
        : override_value;
    }
    else
      value = evaluate(regularization_parameters_.line_distance_threshold, vertex);
    CGAL_precondition(value >= 0 && std::isfinite(value));
    return value;
  }

public:
  Line_distance_threshold_provider(
    const PolygonMesh& pmesh, const MeshPointPropertyMap& mesh_points,
    const ElementAndDomainTranslator& translator,
    const VertexPointMap& vertex_points, const VertexDimensionMap& vertex_dimensions,
    const VertexToMeshElementMap& vertex_to_mesh_element, const Domain& domain,
    const RegularizationParameters& regularization_parameters)
    : pmesh_(pmesh), mesh_points_(mesh_points), translator_(translator),
      vertex_points_(vertex_points), vertex_dimensions_(vertex_dimensions),
      vertex_to_mesh_element_(vertex_to_mesh_element), domain_(domain),
      regularization_parameters_(regularization_parameters)
  {}

  template <typename Vertex>
  double operator()(const Vertex& vertex, Line_status status) const
  {
    if (status == Line_status::isthmus)
      return effective_threshold<
        RegularizationParameters::isthmus_threshold_specified>(
          regularization_parameters_.isthmus_line_distance_threshold, vertex);
    if (status == Line_status::simple)
      return effective_threshold<
        RegularizationParameters::simple_threshold_specified>(
          regularization_parameters_.simple_line_distance_threshold, vertex);
    return effective_threshold<
      RegularizationParameters::corner_threshold_specified>(
        regularization_parameters_.corner_line_distance_threshold, vertex);
  }

  template <typename Vertex>
  double maximum_threshold(const Vertex& vertex) const
  {
    return (std::max)({
      (*this)(vertex, Line_status::isthmus),
      (*this)(vertex, Line_status::simple),
      (*this)(vertex, Line_status::other)
    });
  }
};

} /* namespace internal */

// returns a feature graph extracted from the surface
// PolygonMesh is a model of FaceListGraph
// MeshPointPropertyMap holds the points of the mesh
// ElementAndDomainTranslator is used to:
//   - translate an element descriptor of the mesh to an element descriptor of the sharpness estimator
//   - translate the mesh to the domain of the sharpness estimator
template <typename K, typename PolygonMesh, typename MeshPointPropertyMap,
  typename SharpnessEstimator, typename ElementAndDomainTranslator, typename NormalEstimator,
  typename RegularizationParameters>
internal::Feature_graph<K> detect_sharp_features(
  const PolygonMesh& pmesh, const MeshPointPropertyMap& points,
  const SharpnessEstimator& sharpness_estimator,
  const ElementAndDomainTranslator& element_and_domain_translator,
  const NormalEstimator& normal_estimator,
  const RegularizationParameters& regularization_parameters)
{
  typedef typename K::Point_3 Point_3;
  typedef typename internal::Feature_graph<K> Feature_graph;
  typedef typename internal::Construction_feature_graph Construction_feature_graph;

#ifdef CGAL_FEATURE_GRAPH_DEBUG_OUTPUT_FILES
    internal::debug::output_sharpness_on_surface_elements<K>(
      "surface_sharpness.ply", pmesh, points, sharpness_estimator,
      element_and_domain_translator, normal_estimator);
#endif

  // make initial feature graph containing all sharp elements,
  // with nodes on vertices, edges and faces
  Construction_feature_graph construction_feature_graph;
  auto index_map = get(boost::vertex_index, construction_feature_graph);
  auto vertex_point_map = boost::make_vector_property_map<Point_3>(index_map);
  auto vertex_dimension_map = boost::make_vector_property_map<char>(index_map);
  auto vertex_to_mesh_element_map =
    boost::make_vector_property_map<internal::Face_list_graph_descriptor_holder<PolygonMesh>>(
      index_map);
  auto vertex_sharpness_map =
    boost::make_vector_property_map<typename SharpnessEstimator::Sharpness_number_type>(index_map);
  auto vertex_is_next_to_hole_map = boost::make_vector_property_map<bool>(index_map);
  auto vertex_is_border_map = boost::make_vector_property_map<bool>(index_map);
  make_initial_feature_graph<K>(pmesh, points, sharpness_estimator, element_and_domain_translator,
    construction_feature_graph, vertex_point_map, vertex_dimension_map, vertex_to_mesh_element_map,
    vertex_sharpness_map, vertex_is_next_to_hole_map, vertex_is_border_map);

#ifdef CGAL_FEATURE_GRAPH_DEBUG_OUTPUT_FILES
  {
    internal::debug::output_feature_graph_properties(
      "initial_graph", construction_feature_graph, vertex_point_map, vertex_dimension_map,
      vertex_sharpness_map, vertex_is_next_to_hole_map);
  }
#endif

  // thin the feature graph to lines (pointel isthmus are protected)
  std::vector<bool> vertex_is_isthmus_protected(num_vertices(construction_feature_graph));
  auto vertex_is_isthmus_protected_map =
    boost::make_iterator_property_map(vertex_is_isthmus_protected.begin(),
      get(boost::vertex_index, construction_feature_graph));
  for (auto vertex : boost::make_iterator_range(vertices(construction_feature_graph)))
    put(vertex_is_isthmus_protected_map, vertex, get(vertex_dimension_map, vertex) == 0);
  thin_list_graph_with_priority(
    vertex_sharpness_map, vertex_is_isthmus_protected_map, vertex_is_border_map,
    construction_feature_graph, vertex_is_next_to_hole_map);

  #ifdef CGAL_FEATURE_GRAPH_DEBUG_OUTPUT_FILES
  {
    internal::debug::output_feature_graph_properties(
      "thinned_graph", construction_feature_graph, vertex_point_map, vertex_dimension_map,
      vertex_sharpness_map, vertex_is_next_to_hole_map);
  }
  #endif

  // regularize: remove small or duplicated lines
  const auto regularization_domain = element_and_domain_translator.translate_domain(pmesh);
  internal::Line_distance_threshold_provider<K, PolygonMesh, MeshPointPropertyMap,
    ElementAndDomainTranslator, decltype(vertex_point_map), decltype(vertex_dimension_map),
    decltype(vertex_to_mesh_element_map), decltype(regularization_domain),
    RegularizationParameters> threshold_provider(
      pmesh, points, element_and_domain_translator, vertex_point_map, vertex_dimension_map,
      vertex_to_mesh_element_map, regularization_domain, regularization_parameters);
  internal::Double_line_selector<K, decltype(threshold_provider)> line_selector(threshold_provider);
  regularize_lines(construction_feature_graph, vertex_point_map, line_selector);

#ifdef CGAL_FEATURE_GRAPH_DEBUG_OUTPUT_FILES
  {
    internal::debug::output_feature_graph_properties(
      "regularized_graph", construction_feature_graph, vertex_point_map, vertex_dimension_map,
      vertex_sharpness_map, vertex_is_next_to_hole_map);
  }
#endif

  // [TODO] cut border elements and paste to another graph. Later, they will be merged with other patches

  // optimize [TODO]
  // make output feature_graph [TODO]

  // helper::bBlock("Optimize positions");
  // Polylines<Point> sharp_features_optimized = placement_optimizer(sharp_features_regularized, surface_graph);
  // helper::eBlock();
  // if (!out_filepath.empty())
  // {
  //     Polylines_process::write_polylines<Point, Vector>(sharp_features_optimized,
  //                                                           boost::filesystem::path(out_folder).append("sharp_features_optimized").string());
  // }

  // helper::bBlock("Place polylines in world space");
  // Polylines<Point> sharp_features_world_space = Polylines_process::make_world_polyline<Vector>(sharp_features_optimized, image);
  // helper::eBlock();
  // if (!out_filepath.empty())
  // {
  //     Polylines_process::write_polylines<Point, Vector>(sharp_features_world_space,
  //                                                           boost::filesystem::path(out_folder).append("sharp_features_world_space").string());
  // }

  // return sharp_features_world_space;
  // do nothing for now
  return internal::Feature_graph<K>();
}

// returns a feature graph extracted from the surface
// Surface_graph is a model of FaceListGraph
template <typename K, typename SurfaceGraph, typename CGAL_NP_TEMPLATE_PARAMETERS>
internal::Feature_graph<K> detect_sharp_features(
  const SurfaceGraph& surface, const CGAL_NP_CLASS& np = parameters::default_values())
{
  // do nothing for now
  return internal::Feature_graph<K>();
}

#ifdef CGAL_DGTAL_ENABLED

template <typename K>
struct Detect_sharp_features_on_labeled_image
{
  typedef typename K::Point_3 Point_3;
  typedef std::size_t Size;
  typedef typename K::FT FT;

  typedef internal::Feature_graph<K> Feature_graph;

  template <typename Image, typename CGAL_NP_TEMPLATE_PARAMETERS>
  Feature_graph operator()(
    const Image& image, const CGAL_NP_CLASS& np = parameters::default_values()) const
  {
    typedef Surface_mesh<Point_3> Mesh;
    typedef typename boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;

    using parameters::get_parameter;
    using parameters::get_parameter_reference;
    using parameters::choose_parameter;

    auto value_outside = choose_parameter(get_parameter(np, internal_np::voxel_value), 0);
    // sharpness and normal estimators
    using Sharpness_estimator =
      typename internal_np::Lookup_named_param_def<internal_np::sharpness_estimator_t,
        CGAL_NP_CLASS,
        CGAL::Default>::type;
    using Normal_estimator =
      typename internal_np::Lookup_named_param_def<internal_np::normal_estimator_t,
        CGAL_NP_CLASS,
        CGAL::Default>::type;
    using Isthmus_line_threshold =
      typename internal_np::Lookup_named_param_def<internal_np::isthmus_line_distance_threshold_t,
        CGAL_NP_CLASS, CGAL::Default>::type;
    using Simple_line_threshold =
      typename internal_np::Lookup_named_param_def<internal_np::simple_line_distance_threshold_t,
        CGAL_NP_CLASS, CGAL::Default>::type;
    using Corner_line_threshold =
      typename internal_np::Lookup_named_param_def<internal_np::corner_line_distance_threshold_t,
        CGAL_NP_CLASS, CGAL::Default>::type;
    //  initialize the default estimator only if one (or more) estimator is not given
    Image_AmbrosioTortorelli<K> default_estimator;
    if constexpr (std::is_same<Sharpness_estimator, CGAL::Default>::value ||
                  std::is_same<Normal_estimator, CGAL::Default>::value)
      default_estimator = Image_AmbrosioTortorelli<K>(image);

    auto sharpness_estimator =
      choose_parameter(get_parameter(np, internal_np::sharpness_estimator),
        default_estimator.sharpness_estimator());
    auto normal_estimator =
      choose_parameter(get_parameter(np, internal_np::normal_estimator),
        default_estimator.normal_estimator());
    // regularization parameters
    auto line_distance_threshold =
      choose_parameter(get_parameter(np, internal_np::line_distance_threshold), FT(4.0));
    auto isthmus_line_distance_threshold =
      choose_parameter(get_parameter(np, internal_np::isthmus_line_distance_threshold), FT(-1.0));
    auto simple_line_distance_threshold =
      choose_parameter(get_parameter(np, internal_np::simple_line_distance_threshold), FT(-1.0));
    auto corner_line_distance_threshold =
      choose_parameter(get_parameter(np, internal_np::corner_line_distance_threshold), FT(-1.0));
    constexpr bool has_isthmus_line_threshold =
      !std::is_same<Isthmus_line_threshold, CGAL::Default>::value;
    constexpr bool has_simple_line_threshold =
      !std::is_same<Simple_line_threshold, CGAL::Default>::value;
    constexpr bool has_corner_line_threshold =
      !std::is_same<Corner_line_threshold, CGAL::Default>::value;
    using Regularization_parameters =
      internal::Line_regularization_parameters<has_isthmus_line_threshold,
        has_simple_line_threshold, has_corner_line_threshold,
        decltype(line_distance_threshold), decltype(isthmus_line_distance_threshold),
        decltype(simple_line_distance_threshold), decltype(corner_line_distance_threshold)>;
    const Regularization_parameters regularization_parameters = {
      line_distance_threshold, isthmus_line_distance_threshold,
      simple_line_distance_threshold, corner_line_distance_threshold
    };

    // make a surface graph for each surface of the labeled image
    std::vector<Mesh> meshes;
    std::vector<std::array<float,2>> mesh_domains;
    std::vector<std::map<std::size_t, vertex_descriptor>> image_point_to_mesh_vertex;
    internal::convert_labelled_image_to_face_graphs<K>(
        image, meshes, mesh_domains, image_point_to_mesh_vertex,
        internal::Get_points_property_from_surface_mesh<Point_3>(),
        internal::Get_surface_mesh_domains_and_orientations(),
        value_outside);

    // call detect_sharp_features for each mesh
    internal::Sharpness_estimator_on_image_translator<K, Mesh, Image> sharpness_translator(image);
    std::vector<Feature_graph> sharp_features_for_meshes;
    sharp_features_for_meshes.resize(meshes.size());
    for (std::size_t mesh_index = 0, number_of_mesh = meshes.size();
      mesh_index < number_of_mesh; mesh_index++)
    {
      const Mesh& mesh = meshes[mesh_index];
      sharp_features_for_meshes[mesh_index] =
        detect_sharp_features<K>(mesh, mesh.points(),
          sharpness_estimator, sharpness_translator, normal_estimator,
          regularization_parameters);
    }

    // merge the border lines [TODO]

    // [TODO]do nothing for now
    return Feature_graph();
  }
};

#endif // ifdef CGAL_DGTAL_ENABLED

} /* namespace Feature_graph */

} /* namespace CGAL */

#endif // CGAL_FEATURE_GRAPH_DETECT_SHARP_FEATURES_H