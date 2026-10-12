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

#include <boost/graph/adjacency_list.hpp>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/Surface_mesh/Surface_mesh.h>

#include <CGAL/Feature_graph/internal/feature_graph_operation.h>
#ifdef CGAL_DGTAL_ENABLED
#include <CGAL/Feature_graph/Image_AmbrosioTortorelli.h>
#endif
#include <CGAL/Feature_graph/internal/convert_image_to_surface_mesh.h>

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

  Sharpness_estimator_on_image_translator(const Image& image)
    : image(image)
  { }

  std::size_t get_Khalimsky_coordinate_index(const Point_3& point) const
  {
    double image_space_x = (point[0] - static_cast<double>(image.tx())) / static_cast<double>(image.vx());
    double image_space_y = (point[1] - static_cast<double>(image.ty())) / static_cast<double>(image.vy());
    double image_space_z = (point[2] - static_cast<double>(image.tz())) / static_cast<double>(image.vz());
    // round to nearest coords (+.5) in Khalimsky coordinates (*2)
    int pointel_x = static_cast<int>(2*image_space_x +.5);
    int pointel_y = static_cast<int>(2*image_space_y +.5);
    int pointel_z = static_cast<int>(2*image_space_z +.5);
    // get index
    return pointel_z*image.xdim()*image.ydim()*4 + pointel_y*image.xdim()*2 + pointel_x;
  }

  template <typename PointPropertyMap>
  std::size_t translate_vertex_element(const vertex_descriptor& vertex, const PolygonMesh& /*pmesh*/, const PointPropertyMap& points) const
  {
    return get_Khalimsky_coordinate_index(get(points, vertex));
  }

  template <typename PointPropertyMap>
  std::size_t translate_halfedge_element(const halfedge_descriptor& edge, const PolygonMesh& pmesh, const PointPropertyMap& points) const
  {
    Point_3 midpoint = CGAL::midpoint(get(points, source(edge, pmesh)), get(points, target(edge, pmesh)));
    return get_Khalimsky_coordinate_index(midpoint);
  }

  template <typename PointPropertyMap>
  std::size_t translate_face_element(const face_descriptor& face, const PolygonMesh& pmesh, const PointPropertyMap& points) const
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
    barycenter /= number_of_edges;
    return get_Khalimsky_coordinate_index(CGAL::ORIGIN + barycenter);
  }

  Image translate_domain(const PolygonMesh& /*pmesh*/) const
  {
    return image;
  }
};

} /* namespace internal */

// returns a feature graph extracted from the surface
// PolygonMesh is a model of FaceListGraph
// MeshPointPropertyMap holds the points of the mesh
// ElementAndDomainTranslator is used to:
//   - translate an element descriptor of the mesh to an element descriptor of the sharpness estimator
//   - translate the mesh to the domain of the sharpness estimator
template <typename K,
  typename PolygonMesh, typename MeshPointPropertyMap,
  typename SharpnessEstimator, typename ElementAndDomainTranslator,
  typename NormalEstimator>
internal::Feature_graph<K> detect_sharp_features(
  const PolygonMesh& pmesh, const MeshPointPropertyMap& points,
  const SharpnessEstimator& sharpness_estimator, const ElementAndDomainTranslator& element_and_domain_translator,
  const NormalEstimator& normal_estimator)
{
  typedef typename K::Point_3 Point_3;
  typedef typename K::Vector_3 Vector_3;
  typedef typename internal::Feature_graph<K> Feature_graph;
  typedef typename internal::Construction_feature_graph<K> Construction_feature_graph;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;

  using parameters::get_parameter;
  using parameters::get_parameter_reference;
  using parameters::choose_parameter;

  // make initial feature graph
  Construction_feature_graph construction_feature_graph;
  boost::vector_property_map<char> feature_graph_element_dimension;
  boost::vector_property_map<internal::Face_list_graph_descriptor_holder<PolygonMesh>> feature_graph_element_to_mesh_element_property_map;
  boost::vector_property_map<typename SharpnessEstimator::Sharpness_number_type> feature_graph_sharpness_property_map;
  boost::vector_property_map<bool> feature_graph_next_to_hole_property_map;
  make_initial_feature_graph<K>(pmesh, points, sharpness_estimator, element_and_domain_translator,
    construction_feature_graph, feature_graph_element_dimension, feature_graph_element_to_mesh_element_property_map,
    feature_graph_sharpness_property_map, feature_graph_next_to_hole_property_map);

  // {
  //   std::ofstream points_output("points.xyz");
  //   std::ofstream edges_output("edges.polylines.txt");
  //   for (auto edge_it = edges(feature_graph).first, edge_end = edges(feature_graph).second; edge_it != edge_end; ++edge_it)
  //   {
  //     auto edge = *edge_it;
  //     auto vertex = source(edge, feature_graph);
  //     int dim = static_cast<float>(get(feature_graph_element_dimension, vertex)) * 255.0 / 3.0;
  //     int sharpness = static_cast<float>(get(feature_graph_sharpness_property_map, vertex)) * 255.0;
  //     points_output << feature_graph[vertex] << " " << dim << " " << sharpness << " " << sharpness << std::endl;

  //     edges_output << "2 " << feature_graph[vertex] << " " << feature_graph[target(edge, feature_graph)] << std::endl;
  //   }
  //   points_output.close();
  //   edges_output.close();
  // }

  // thin to lines

  // regularize
  // optimize

  // // Make initial feature graph
  // Polylines<Point> sharp_features = polyline_extractor.template operator()<Point, Vector>(sharp_surface_graph, deleted_faces, pointel_values);
  // helper::eBlock();
  // if (!out_filepath.empty())
  // {
  //     Polylines_process::write_polylines<Point, Vector>(sharp_features,
  //                                                           boost::filesystem::path(out_folder).append("sharp_features_initial").string());
  // }

  // helper::bBlock("Regularize double lines");
  // Polylines<Point> sharp_features_regularized = regularization_functor(sharp_features);
  // helper::eBlock();
  // if (!out_filepath.empty())
  // {
  //     Polylines_process::write_polylines<Point, Vector>(sharp_features_regularized,
  //                                                           boost::filesystem::path(out_folder).append("sharp_features_regularized").string());
  // }

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
internal::Feature_graph<K> detect_sharp_features(const SurfaceGraph& surface, const CGAL_NP_CLASS& np = parameters::default_values())
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
  Feature_graph operator()(const Image& image, const CGAL_NP_CLASS& np = parameters::default_values()) const
  {
    typedef Surface_mesh<Point_3> Mesh;
    typedef typename boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;

    using parameters::get_parameter;
    using parameters::get_parameter_reference;
    using parameters::choose_parameter;

    auto value_outside = choose_parameter(get_parameter(np, internal_np::voxel_value), 0);
    // sharpness and normal estimators
    using Sharpness_estimator = typename internal_np::Lookup_named_param_def<internal_np::sharpness_estimator_t,
                                                                                  CGAL_NP_CLASS,
                                                                                  CGAL::Default>::type;
    using Normal_estimator = typename internal_np::Lookup_named_param_def<internal_np::normal_estimator_t,
                                                                                  CGAL_NP_CLASS,
                                                                                  CGAL::Default>::type;
    //  initialize the default estimator only if one (or more) estimator is not given
    Image_AmbrosioTortorelli<K> default_estimator;
    if constexpr (std::is_same<Sharpness_estimator, CGAL::Default>::value ||
                  std::is_same<Normal_estimator, CGAL::Default>::value)
      default_estimator = Image_AmbrosioTortorelli<K>(image);

    auto sharpness_estimator = choose_parameter(get_parameter(np, internal_np::sharpness_estimator), default_estimator.sharpness_estimator());
    auto normal_estimator = choose_parameter(get_parameter(np, internal_np::normal_estimator), default_estimator.normal_estimator());

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
    for (std::size_t mesh_index = 0, number_of_mesh = meshes.size(); mesh_index < number_of_mesh; mesh_index++)
    {
      const Mesh& mesh = meshes[mesh_index];
      sharp_features_for_meshes[mesh_index] = detect_sharp_features<K>(mesh, mesh.points(), sharpness_estimator, sharpness_translator, normal_estimator);
    }

    // merge the border lines

    // do nothing for now
    return Feature_graph();
  }
};

#endif // ifdef CGAL_DGTAL_ENABLED

} /* namespace Feature_graph */

} /* namespace CGAL */

#endif // CGAL_FEATURE_GRAPH_DETECT_SHARP_FEATURES_H