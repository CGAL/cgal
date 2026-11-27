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

#ifndef CGAL_FEATURE_GRAPH_INTERNAL_FEATURE_GRAPH_OPERATION_H
#define CGAL_FEATURE_GRAPH_INTERNAL_FEATURE_GRAPH_OPERATION_H

#include <CGAL/license/Feature_graph.h>

#include <cstddef>
#include <map>
#include <variant>
#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>

#include <CGAL/boost/graph/helpers.h>

namespace CGAL
{

namespace Feature_graph
{

namespace internal
{

// feature graph type returned by detect_sharp_features. Optimized for fast iterations.
template <typename K>
using Feature_graph =
  boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS, typename K::Point_3>;

// feature graph type used for at the creation.
// Optimized for fast addition and removal of elements, but do not store points.
// Preserve vertex and edge iterators on addition and removal
using Construction_feature_graph =
  boost::adjacency_list<boost::listS, boost::listS, boost::undirectedS,
    boost::property<boost::vertex_index_t, std::size_t>>;

template <typename FaceListGraph>
class Face_list_graph_descriptor_holder : public std::variant<
  typename boost::graph_traits<FaceListGraph>::vertex_descriptor,
  typename boost::graph_traits<FaceListGraph>::halfedge_descriptor,
  typename boost::graph_traits<FaceListGraph>::face_descriptor>
{
  typedef typename std::variant<
    typename boost::graph_traits<FaceListGraph>::vertex_descriptor,
    typename boost::graph_traits<FaceListGraph>::halfedge_descriptor,
    typename boost::graph_traits<FaceListGraph>::face_descriptor> Base;
  using Base::variant;
};

} /* namespace internal */

// returns a feature graph containing the selected elements
//
//  Selected PolygonMesh facet     Feature_graph equivalent
//          X <---- X                     X - X - X
//          |       ^                     | \ | / |
//         ...      |         -->        ...- X - X
//          V       |                     | / | \ |
//          X -fe-> X                     X - X - X
//
// PolygonMesh is model of FaceListGraph
// ElementAndDomainTranslator is used to:
//   - translate an element descriptor of the mesh to an element descriptor of the sharpness estimator
//   - translate the mesh to the domain of the sharpness estimator
// VertexPointMap assigns a `K::Point_3` point to each feature graph element.
// VertexDimensionPropertyMap assigns a `DimensionTag` of either 0, 1 or 2 to each feature graph element.
// VertexToMeshElementPropertyMap assigns a `Face_list_graph_descriptor_holder<PolygonMesh>`
//   to each feature graph element.
// VertexSharpnessPropertyMap assigns a `SharpnessEstimator::Sharpness_value_type` to each feature graph element
// VertexNextToHolePropertyMap assigns true to each feature graph element that is next to a hole
//   (i.e. a face that is not sharp).
// VertexBorderPropertyMap assigns true to each feature graph element that is border element of the surface.
template <typename K,
  typename PolygonMesh, typename MeshPointPropertyMap,
  typename SharpnessEstimator, typename ElementAndDomainTranslator,
  typename VertexPointMap, typename VertexDimensionPropertyMap,
  typename VertexToMeshElementPropertyMap, typename VertexSharpnessPropertyMap,
  typename VertexNextToHolePropertyMap, typename VertexBorderPropertyMap>
void make_initial_feature_graph(
  const PolygonMesh& pmesh, const MeshPointPropertyMap& points,
  const SharpnessEstimator& sharpness_estimator,
  const ElementAndDomainTranslator& element_and_domain_translator,
  internal::Construction_feature_graph& feature_graph,
  VertexPointMap& vertex_point_map,
  VertexDimensionPropertyMap& vertex_dimension_map,
  VertexToMeshElementPropertyMap& vertex_to_mesh_element_map,
  VertexSharpnessPropertyMap& vertex_sharpness_map,
  VertexNextToHolePropertyMap& vertex_is_next_to_hole_map,
  VertexBorderPropertyMap& vertex_is_border_map)
{
  typedef typename K::Point_3 Point_3;
  typedef typename K::Vector_3 Vector_3;
  typedef internal::Construction_feature_graph Feature_graph;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor edge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
  typedef typename internal::Face_list_graph_descriptor_holder<PolygonMesh> Face_list_graph_descriptor_holder;

  typedef typename boost::graph_traits<Feature_graph>::vertex_descriptor Feature_graph_element;

  auto sharpness_estimator_domain = element_and_domain_translator.translate_domain(pmesh);

  std::map<vertex_descriptor, Feature_graph_element> vertex_to_element;
  std::map<edge_descriptor, Feature_graph_element> edge_to_element;
  std::vector<face_descriptor> smooth_faces;
  // for each face
  for (const face_descriptor& face : faces(pmesh))
  {
    // get or create elements for any polygon
    /*  X <---- X         X - X - X
    *   |       ^         | \ | / |
    *  ...      |   -->  ...- X - X
    *   V       |         | / | \ |
    *   X -fe-> X         X - X - X
    */
    halfedge_descriptor face_edge = halfedge(face, pmesh);
    // filter non sharp faces
    auto face_element_sharpness =
      sharpness_estimator.template operator()<Dimension_tag<2>>(
        element_and_domain_translator.translate_face_element(face, pmesh, points),
        sharpness_estimator_domain);
    if (face_element_sharpness < 0)
    {
      smooth_faces.push_back(face);
      continue;
    }
    Feature_graph_element face_element = add_vertex(num_vertices(feature_graph), feature_graph);
    put(vertex_dimension_map, face_element, 2);
    Face_list_graph_descriptor_holder face_variant = face;
    put(vertex_to_mesh_element_map, face_element, face_variant);
    put(vertex_sharpness_map, face_element, face_element_sharpness);
    put(vertex_is_next_to_hole_map, face_element, false);
    Vector_3 barycenter(0,0,0);
    int number_of_edges = 0;
    Feature_graph_element previous_graph_edge;
    bool previous_graph_edge_already_created;
    Feature_graph_element first_feature_graph_source_vertex;
    bool first_source_vertex_already_created;
    // iterate on halfedge around the face (circulator)
    halfedge_descriptor current_face_edge = face_edge;
    do
    {
      // get or create source vertex
      vertex_descriptor source_vertex = source(current_face_edge, pmesh);
      Feature_graph_element feature_graph_source_vertex;
      auto source_vertex_element_it = vertex_to_element.find(source_vertex);
      bool is_source_vertex_already_created = source_vertex_element_it != vertex_to_element.end();
      if (!is_source_vertex_already_created)
      {
        feature_graph_source_vertex = add_vertex(num_vertices(feature_graph), feature_graph);
        vertex_to_element[source_vertex] = feature_graph_source_vertex;
        put(vertex_point_map, feature_graph_source_vertex, get(points, source_vertex));
        put(vertex_dimension_map, feature_graph_source_vertex, 0);
        put(vertex_to_mesh_element_map, feature_graph_source_vertex,
          Face_list_graph_descriptor_holder(source_vertex));
        auto source_vertex_sharpness =
          sharpness_estimator.template operator()<Dimension_tag<0>>(
            element_and_domain_translator.translate_vertex_element(source_vertex, pmesh, points),
            sharpness_estimator_domain);
        put(vertex_sharpness_map, feature_graph_source_vertex, source_vertex_sharpness);
        put(vertex_is_next_to_hole_map, feature_graph_source_vertex, false);
        put(vertex_is_border_map, feature_graph_source_vertex, false);
      }
      else
        feature_graph_source_vertex = source_vertex_element_it->second;
      if (number_of_edges == 0)
      {
        first_feature_graph_source_vertex = feature_graph_source_vertex;
        first_source_vertex_already_created = is_source_vertex_already_created;
      }
      // get or create edge (use edge to avoid duplicates)
      edge_descriptor current_as_edge = edge(current_face_edge, pmesh);
      Feature_graph_element feature_graph_edge;
      auto current_edge_element_it = edge_to_element.find(current_as_edge);
      bool is_edge_already_created = current_edge_element_it != edge_to_element.end();
      if (!is_edge_already_created)
      {
        feature_graph_edge = add_vertex(num_vertices(feature_graph), feature_graph);
        edge_to_element[current_as_edge] = feature_graph_edge;
        Point_3 edge_point =
          CGAL::midpoint(get(points, source_vertex), get(points, target(current_face_edge, pmesh)));
        put(vertex_point_map, feature_graph_edge, edge_point);
        put(vertex_dimension_map, feature_graph_edge, 1);
        put(vertex_to_mesh_element_map, feature_graph_edge,
          Face_list_graph_descriptor_holder(current_face_edge));
        auto edge_sharpness =
          sharpness_estimator.template operator()<Dimension_tag<1>>(
            element_and_domain_translator.translate_halfedge_element(current_face_edge, pmesh, points),
            sharpness_estimator_domain);
        put(vertex_sharpness_map, feature_graph_edge, edge_sharpness);
        put(vertex_is_next_to_hole_map, feature_graph_edge, false);
        put(vertex_is_border_map, feature_graph_edge, false);
      }
      else
        feature_graph_edge = current_edge_element_it->second;
      if (CGAL::is_border_edge(current_face_edge, pmesh))
      {
        put(vertex_is_border_map, feature_graph_source_vertex, true);
        put(vertex_is_border_map, feature_graph_edge, true);
      }
      // make links
      add_edge(face_element, feature_graph_edge, feature_graph);
      add_edge(face_element, feature_graph_source_vertex, feature_graph);
      if (!is_source_vertex_already_created || !is_edge_already_created)
        add_edge(feature_graph_source_vertex, feature_graph_edge, feature_graph);
      if (number_of_edges > 0 &&
        (!is_source_vertex_already_created || !previous_graph_edge_already_created))
        add_edge(feature_graph_source_vertex, previous_graph_edge, feature_graph);
      // iterate
      barycenter += Vector_3(CGAL::ORIGIN, get(points, source_vertex));
      number_of_edges++;
      current_face_edge = next(current_face_edge, pmesh);
      previous_graph_edge = feature_graph_edge;
      previous_graph_edge_already_created = is_edge_already_created;
    } while (current_face_edge != face_edge);
    if (number_of_edges != 0)
      barycenter /= number_of_edges;
    put(vertex_point_map, face_element, CGAL::ORIGIN + barycenter);
    if (number_of_edges > 0 &&
      (!first_source_vertex_already_created || !previous_graph_edge_already_created))
      add_edge(first_feature_graph_source_vertex, previous_graph_edge, feature_graph);
  }
  // update "next_to_hole" property of elements for all feature graph element adjacent to smooth faces
  for (const face_descriptor& face : smooth_faces)
  {
    halfedge_descriptor face_edge = halfedge(face, pmesh);
    // iterate on halfedge around the face
    halfedge_descriptor current_face_edge = face_edge;
    do
    {
      vertex_descriptor source_vertex = source(current_face_edge, pmesh);
      auto feature_graph_source_vertex_iterator = vertex_to_element.find(source_vertex);
      if (feature_graph_source_vertex_iterator != vertex_to_element.end())
        put(vertex_is_next_to_hole_map, feature_graph_source_vertex_iterator->second, true);
      edge_descriptor current_as_edge = edge(current_face_edge, pmesh);
      auto feature_graph_edge_iterator = edge_to_element.find(current_as_edge);
      if (feature_graph_edge_iterator != edge_to_element.end())
        put(vertex_is_next_to_hole_map, feature_graph_edge_iterator->second, true);
      current_face_edge = next(current_face_edge, pmesh);
    } while (current_face_edge != face_edge);
  }
}

} /* namespace Feature_graph */

} /* namespace CGAL */

#endif // CGAL_FEATURE_GRAPH_INTERNAL_FEATURE_GRAPH_OPERATION_H