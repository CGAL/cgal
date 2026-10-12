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

#ifndef CGAL_FEATURE_GRAPH_INTERNAL_FEATURE_GRAPH_H
#define CGAL_FEATURE_GRAPH_INTERNAL_FEATURE_GRAPH_H

#include <CGAL/license/Feature_graph.h>

#include <boost/graph/adjacency_list.hpp>

namespace CGAL
{

namespace Feature_graph
{

namespace internal
{

// The feature graph returned by detect_sharp_features. Optimized for fast iterations.
template <typename K>
class Feature_graph : public boost::adjacency_list<boost::listS,boost::vecS,boost::undirectedS,typename K::Point_3>
{
  typedef typename boost::adjacency_list<boost::listS,boost::vecS,boost::undirectedS,typename K::Point_3> Base;
  using Base::adjacency_list;
};

// The feature graph used for at the creation. Optimized for fast addition and removals.
template <typename K>
class Construction_feature_graph : public boost::adjacency_list<boost::listS,boost::listS,boost::undirectedS,typename K::Point_3>
{
  typedef typename boost::adjacency_list<boost::listS,boost::listS,boost::undirectedS,typename K::Point_3> Base;
  using Base::adjacency_list;
};

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

// delete the first removable element from the graph
// start search from "lowest_priority_to_check_iterator"
template <typename K, typename PriorityValue, typename IthmusProtectedPropertyMap, typename NextToHolePropertyMap, typename Graph>
void thin_remove_lowest_value_vertex(
  const PriorityValue& priority_value, const IthmusProtectedPropertyMap& ithmus_protected_property_map,
  internal::Construction_feature_graph<K>& feature_graph,
  NextToHolePropertyMap& next_to_hole_property_map,
  std::map<double, std::list<internal::Construction_feature_graph<K>::vertex_descriptor>>& sorted_elements,
  std::map<double, std::list<internal::Construction_feature_graph<K>::vertex_descriptor>>::iterator& lowest_priority_to_check_iterator)
{
  typedef internal::Construction_feature_graph<K> Feature_graph;
  using Vertex = Graph::vertex_descriptor;

  using Element_graph = Element_graph<Point>;
  using Graph_element = typename Element_graph::Graph_element;
  using Element_iterator = typename Element_graph::Element_iterator;
  using Prioriry_iterator = typename std::map<double, std::list<Element_iterator>>::iterator;

  typename std::list<Element_iterator>::iterator value_element_it;
  typename std::list<Element_iterator>::iterator value_element_end;

  for (Prioriry_iterator sorted_element_end = priority_map.end(); priority_begin != sorted_element_end; priority_begin++)
  {
      std::list<Element_iterator>& value_elements = priority_begin->second;

      for (value_element_it = value_elements.begin(), value_element_end = value_elements.end();
          value_element_it != value_element_end; ++value_element_it)
      {
          Element_iterator& element_it = *value_element_it;
          Graph_element* element = element_it->get();
          if (is_element_removable(graph, element))
          {
              // update the start of search (priority_begin) to the lowest neighbor
              const double& element_value = priority(graph, *element);
              double lowest_neighbor_value = element_value;
              for ( auto& element_neigh : element->neighbors )
              {
                  const double& element_value = priority(graph, *element_neigh);
                  if (element_value < lowest_neighbor_value)
                  {
                      lowest_neighbor_value = element_value;
                  }
              }
              // erase from graph
              graph.remove_element(element_it);
              // erase from priority map
              value_elements.erase(value_element_it);
              // if (value_elements.size() == 0)
              //     priority_map.erase(priority_begin);
              if (lowest_neighbor_value != element_value)
              {
                  priority_begin = priority_map.find(lowest_neighbor_value);
              }
              return true;
          }
      }
  }

  return false;
}

} /* namespace internal */

// returns a feature graph containing the selected elements
// PolygonMesh is model of FaceListGraph
// ElementAndDomainTranslator is used to:
//   - translate an element descriptor of the mesh to an element descriptor of the sharpness estimator
//   - translate the mesh to the domain of the sharpness estimator
// FeatureGraphElementDimensionPropertyMap assignes a dimension of either 0, 1 or 2 to each element of the feature graph.
//      It must define `put(FeatureGraphElementDimensionPropertyMap map, Feature_graph::vertex_descriptor element, DimensionTag dimension)`
// FeatureGraphElementToMeshElementPropertyMap assigns a `Face_list_graph_descriptor_holder` to the polygon mesh for each element of the feature graph.
//      It must define `put(FeatureGraphElementToMeshElementPropertyMap map, Feature_graph::vertex_descriptor element, Face_list_graph_descriptor_holder descriptor)`
// FeatureGraphSharpnessPropertyMap assigns a `SharpnessEstimator::Sharpness_value_type` for each element of the feature graph.
//      It must define `put(FeatureGraphSharpnessPropertyMap map, Feature_graph::vertex_descriptor element, SharpnessEstimator::Sharpness_value_type sharpness)`
// FeatureGraphNextToHolePropertyMap assigns true for each element of the feature graph that are next to a hole (i.e. a face that is not sharp).
//      It must define `put(FeatureGraphNextToHolePropertyMap map, Feature_graph::vertex_descriptor element, bool is_next_to_hole)`
template <typename K,
  typename PolygonMesh, typename MeshPointPropertyMap,
  typename SharpnessEstimator, typename ElementAndDomainTranslator,
  typename FeatureGraphElementDimensionPropertyMap, typename FeatureGraphElementToMeshElementPropertyMap,
  typename FeatureGraphSharpnessPropertyMap, typename FeatureGraphNextToHolePropertyMap>
void make_initial_feature_graph(
  const PolygonMesh& pmesh, const MeshPointPropertyMap& points,
  const SharpnessEstimator& sharpness_estimator, const ElementAndDomainTranslator& element_and_domain_translator,
  internal::Construction_feature_graph<K>& feature_graph,
  FeatureGraphElementDimensionPropertyMap& feature_graph_element_dimension,
  FeatureGraphElementToMeshElementPropertyMap& feature_graph_element_to_mesh_element_property_map,
  FeatureGraphSharpnessPropertyMap& feature_graph_sharpness_property_map,
  FeatureGraphNextToHolePropertyMap& feature_graph_next_to_hole_property_map)
{
  typedef typename K::Point_3 Point_3;
  typedef typename K::Vector_3 Vector_3;
  typedef internal::Construction_feature_graph<K> Feature_graph;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor edge_descriptor; // TODO avoid duplicate edges
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
  typedef typename internal::Face_list_graph_descriptor_holder<PolygonMesh> Face_list_graph_descriptor_holder;

  typedef typename Feature_graph::vertex_descriptor Freature_graph_element;

  auto sharpness_estimator_domain = element_and_domain_translator.translate_domain(pmesh);

  std::map<vertex_descriptor, Freature_graph_element> vertex_to_element;
  std::map<edge_descriptor, Freature_graph_element> edge_to_element;
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
    auto face_element_sharpness = sharpness_estimator.template operator()<Dimension_tag<2>>(element_and_domain_translator.translate_face_element(face, pmesh, points), sharpness_estimator_domain);
    if (face_element_sharpness < 0)
    {
      smooth_faces.push_back(face);
      continue;
    }
    Freature_graph_element face_element = add_vertex(feature_graph);
    put(feature_graph_element_dimension, face_element, 2);
    Face_list_graph_descriptor_holder face_variant = face;
    put(feature_graph_element_to_mesh_element_property_map, face_element, face_variant);
    put(feature_graph_sharpness_property_map, face_element, face_element_sharpness);
    put(feature_graph_next_to_hole_property_map, face_element, false);
    Vector_3 barycenter(0,0,0);
    int number_of_edges = 0;
    // iterate on halfedge around the face
    halfedge_descriptor current_face_edge = face_edge;
    do
    {
      // get or create source vertex
      vertex_descriptor source_vertex = source(current_face_edge, pmesh);
      Freature_graph_element feature_graph_source_vertex;
      bool is_source_vertex_created = vertex_to_element.find(source_vertex) == vertex_to_element.end();
      if (is_source_vertex_created)
      {
        feature_graph_source_vertex = add_vertex(feature_graph);
        feature_graph[feature_graph_source_vertex] = get(points, source_vertex);
        put(feature_graph_element_dimension, feature_graph_source_vertex, 0);
        put(feature_graph_element_to_mesh_element_property_map, feature_graph_source_vertex, Face_list_graph_descriptor_holder(source_vertex));
        auto source_vertex_sharpness = sharpness_estimator.template operator()<Dimension_tag<0>>(element_and_domain_translator.translate_vertex_element(source_vertex, pmesh, points), sharpness_estimator_domain);
        put(feature_graph_sharpness_property_map, feature_graph_source_vertex, source_vertex_sharpness);
        put(feature_graph_next_to_hole_property_map, feature_graph_source_vertex, false);
      }
      // get or create edge (use edge to avoid duplicates)
      edge_descriptor current_as_edge = edge(current_face_edge, pmesh);
      Freature_graph_element feature_graph_edge;
      bool is_edge_created = edge_to_element.find(current_as_edge) == edge_to_element.end();
      if (is_edge_created)
      {
        feature_graph_edge = add_vertex(feature_graph);
        feature_graph[feature_graph_edge] = CGAL::midpoint(get(points, source_vertex), get(points, target(current_as_edge, pmesh)));
        put(feature_graph_element_dimension, feature_graph_edge, 1);
        put(feature_graph_element_to_mesh_element_property_map, feature_graph_edge, Face_list_graph_descriptor_holder(face_edge));
        auto edge_sharpness = sharpness_estimator.template operator()<Dimension_tag<1>>(element_and_domain_translator.translate_halfedge_element(face_edge, pmesh, points), sharpness_estimator_domain);
        put(feature_graph_sharpness_property_map, feature_graph_edge, edge_sharpness);
        put(feature_graph_next_to_hole_property_map, feature_graph_edge, false);
      }
      // make links
      add_edge(face_element, feature_graph_edge, feature_graph);
      add_edge(face_element, feature_graph_source_vertex, feature_graph);
      if (is_source_vertex_created || is_edge_created)
        add_edge(feature_graph_source_vertex, feature_graph_edge, feature_graph);
      // iterate
      barycenter += Vector_3(CGAL::ORIGIN, get(points, source_vertex));
      number_of_edges++;
      current_face_edge = next(current_face_edge, pmesh);
    } while (current_face_edge != face_edge);
    barycenter /= number_of_edges;
    feature_graph[face_element] = CGAL::ORIGIN + barycenter;
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
        put(feature_graph_next_to_hole_property_map, feature_graph_source_vertex_iterator->second, true);
      edge_descriptor current_as_edge = edge(current_face_edge, pmesh);
      auto feature_graph_edge_iterator = edge_to_element.find(current_as_edge);
      if (feature_graph_edge_iterator != edge_to_element.end())
        put(feature_graph_next_to_hole_property_map, feature_graph_edge_iterator->second, true);
      current_face_edge = next(current_face_edge, pmesh);
    } while (current_face_edge != face_edge);
  }
}

// thins the graph using a priority based on values: lower values are deleted first.
// PriorityValue returns a priority value for each graph element
// IthmusProtectedPropertyMap returns true iff the graph element should be protected when it is isthmus
// NextToHolePropertyMap assigns true if the element is next to a deleted element (i.e. a )
template <typename K, typename PriorityValue, typename IthmusProtectedPropertyMap, typename NextToHolePropertyMap>
void thin_list_graph_with_priority(
  const PriorityValue& priority_value, const IthmusProtectedPropertyMap& ithmus_protected_property_map,
  internal::Construction_feature_graph<K>& feature_graph,
  NextToHolePropertyMap& next_to_hole_property_map)
{
  typedef internal::Construction_feature_graph<K> Feature_graph;
  using Vertex = Graph::vertex_descriptor;
  bool deleted_simple = false;

  // sort vertices by priority values, and create lists for elements with the same values.
  std::map<double, std::list<Vertex>> sorted_elements;
  for (Vertex vertex : vertices(feature_graph))
  {
    sorted_elements[priority_value(vertex, feature_graph)].push_back(vertex);
  }
  typename std::map<double, std::list<Vertex>>::iterator lowest_priority_to_check_iterator = sorted_elements.begin();

  // mark_in_h(neigh(v)); clear_vertex(v, g); remove_vertex(v, g); // does not invalidates desciptor for Construction_feature_graph.

  do
  {
    deleted_simple = internal::thin_remove_lowest_value_vertex( feature_graph, sorted_elements, lowest_priority_to_check_iterator, priority );
  }
  while (deleted_simple);
}

} /* namespace Feature_graph */

} /* namespace CGAL */

#endif // CGAL_FEATURE_GRAPH_INTERNAL_FEATURE_GRAPH_H