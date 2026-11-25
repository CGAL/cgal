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


#ifndef CGAL_FEATURE_GRAPH_INTERNAL_GRAPH_OPERATION_H
#define CGAL_FEATURE_GRAPH_INTERNAL_GRAPH_OPERATION_H

#include <cstddef>
#include <list>
#include <map>
#include <set>
#include <vector>

#include <CGAL/license/Feature_graph.h>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/property_map/vector_property_map.hpp>
#include <boost/range/iterator_range.hpp>

#include <CGAL/Feature_graph/internal/feature_graph_operation.h>


namespace CGAL
{

namespace Feature_graph
{

namespace internal
{

struct CountDiscoveredVerticesVisitor : public boost::default_bfs_visitor {
  std::size_t& number_of_vertices_ref;
  CountDiscoveredVerticesVisitor(std::size_t& number_of_vertices_to_change)
    : number_of_vertices_ref(number_of_vertices_to_change)
  {}

  template <typename Graph>
  void discover_vertex(typename Graph::vertex_descriptor /*v*/, const Graph& /*g*/) {
      number_of_vertices_ref++;
  }
};

// assumes undirected graph with no self edges (v->v not possible)
template <typename Graph>
bool is_vertex_locally_connected(typename Graph::vertex_descriptor vertex, const Graph& graph)
{
  using Vertex = typename Graph::vertex_descriptor;
  // create neighborhood graph
  std::set<Vertex> neighbor_set;
  for (auto out_edge : boost::make_iterator_range(out_edges(vertex, graph)))
    neighbor_set.insert(target(out_edge, graph));
  std::size_t neighbor_size = neighbor_set.size();
  if (neighbor_size <= 1)
    return true;
  auto is_neighbor_predicate = boost::is_in_subset<std::set<Vertex>>(neighbor_set);
  auto induced_subgraph = boost::make_filtered_graph(graph, boost::keep_all(), is_neighbor_predicate);
  // check if it is connected
  std::map<Vertex, boost::default_color_type> colors;
  auto color_map = make_assoc_property_map(colors);
  std::size_t number_of_visited_vertices = 0;
  CountDiscoveredVerticesVisitor vertex_count_visitor(number_of_visited_vertices);
  boost::breadth_first_search(
      induced_subgraph,
      *neighbor_set.begin(),
      boost::visitor(vertex_count_visitor)
            .color_map(color_map)
  );
  return number_of_visited_vertices == neighbor_size;
  // boost::filtered_graph<Graph, boost::keep_all, boost::is_in_subset<std::unordered_set<Vertex>>> filter_graph(graph, boost::keep_all{}, vertex_predicate);
}

// returns true iif the vertex is removable, i.e. it is:
//  - not protected : more than 1 neighbor or not a pointel
//  - is simple :
//      - do not create hole : is next to a hole (a neighbor element has been deleted)
//      - do not merge hole : the subgraph of neighbors is connected
template <typename IsthmusProtectedPropertyMap, typename NextToHolePropertyMap, typename ProtectedPropertyMap>
bool is_vertex_removable(
  const Construction_feature_graph& graph, // [TODO] template
  typename Construction_feature_graph::vertex_descriptor vertex,
  const IsthmusProtectedPropertyMap& isthmus_protected_property_map,
  const NextToHolePropertyMap& next_to_hole_property_map,
  const ProtectedPropertyMap& protected_vertices)
{
  typedef Construction_feature_graph Graph;
  using Vertex = typename Graph::vertex_descriptor;
  // check if the vertex is protected
  if (get(protected_vertices, vertex))
    return false;
  // check if the vertex is isthmus-protected
  if (get(isthmus_protected_property_map, vertex) && out_degree(vertex, graph) == 1)
    return false;
  // check if a hole is created
  if (!get(next_to_hole_property_map, vertex))
    return false;
  // check if holes are merged
  if (!is_vertex_locally_connected(vertex, graph))
    return false;
  // all checks passed
  return true;
}

// delete the first removable element from the graph
// start search from "lowest_priority_to_check_iterator"
template <typename PriorityValuePropertyMap, typename IsthmusProtectedPropertyMap, typename ProtectedPropertyMap,
          typename NextToHolePropertyMap, typename IsComputedPropertyMap>
bool thin_remove_lowest_value_vertex(
  const PriorityValuePropertyMap& priority_value_property_map, const IsthmusProtectedPropertyMap& isthmus_protected_property_map,
  const ProtectedPropertyMap& protected_vertices,
  Construction_feature_graph& graph, // [TODO] template
  NextToHolePropertyMap& next_to_hole_property_map,
  std::map<typename boost::property_traits<PriorityValuePropertyMap>::value_type,
      std::list<typename Construction_feature_graph::vertex_descriptor>>& sorted_elements, // [TODO] replace with sorted list. Requires recursive implementation for efficient checks
  typename std::map<typename boost::property_traits<PriorityValuePropertyMap>::value_type,
      std::list<typename Construction_feature_graph::vertex_descriptor>>::iterator& lowest_priority_to_check_iterator,
  IsComputedPropertyMap& is_element_computed_cache)
{
  typedef Construction_feature_graph Graph;
  using Vertex = typename Graph::vertex_descriptor;
  using Sharpness_type = typename boost::property_traits<PriorityValuePropertyMap>::value_type;

  using Prioriry_iterator = typename std::map<Sharpness_type, std::list<Vertex>>::iterator;

  typename std::list<Vertex>::iterator value_element_it;
  typename std::list<Vertex>::iterator value_element_end;

  for (Prioriry_iterator sorted_element_end = sorted_elements.end(); lowest_priority_to_check_iterator != sorted_element_end; lowest_priority_to_check_iterator++)
  {
    std::list<Vertex>& value_elements = lowest_priority_to_check_iterator->second;

    for (value_element_it = value_elements.begin(), value_element_end = value_elements.end();
      value_element_it != value_element_end; ++value_element_it)
    {
      Vertex vertex = *value_element_it;
      // cheack if already computed that the vertex is not removable
      if (get(is_element_computed_cache, vertex))
        continue;
      if (is_vertex_removable(graph, vertex, isthmus_protected_property_map, next_to_hole_property_map, protected_vertices))
      {
        // update the start of search (lowest_priority_to_check_iterator) to the neighbor with lowest priority
        // also set next_to_hole_property_map to true for all neighbors
        // also remove the cached removable property for all neighbors
        const Sharpness_type& current_vertex_priority = lowest_priority_to_check_iterator->first;
        Sharpness_type lowest_neighbor_priority = current_vertex_priority;
        for (auto out_edge : boost::make_iterator_range(out_edges(vertex, graph)))
        {
          Vertex neighbor_vertex = target(out_edge, graph);
          put(next_to_hole_property_map, neighbor_vertex, true);
          put(is_element_computed_cache, neighbor_vertex, false);
          const Sharpness_type& neighbor_priority = get(priority_value_property_map, neighbor_vertex);
          if (neighbor_priority < lowest_neighbor_priority)
              lowest_neighbor_priority = neighbor_priority;
        }
        // erase from graph
        clear_vertex(vertex, graph);
        remove_vertex(vertex, graph);
        // erase from priority map
        value_elements.erase(value_element_it);
        if (value_elements.empty())
          lowest_priority_to_check_iterator = sorted_elements.erase(lowest_priority_to_check_iterator);
        if (lowest_neighbor_priority < current_vertex_priority)
          lowest_priority_to_check_iterator = sorted_elements.find(lowest_neighbor_priority);
        return true;
      }
      else
        put(is_element_computed_cache, vertex, true);
    }
  }
  // there are no vertex left to remove
  return false;
}


} /* namespace internal */

// thins the graph using a priority based on values: lower values are deleted first.
// PriorityValuePropertyMap returns a priority value for each graph element
// IsthmusProtectedPropertyMap returns true iff the graph element should be protected when it is isthmus
// NextToHolePropertyMap assigns true if the element is next to a deleted element (i.e. a )
// [TODO] make unit test to prove that it works because it is hard to digest.
template <typename PriorityValuePropertyMap, typename IsthmusProtectedPropertyMap, typename ProtectedPropertyMap, typename NextToHolePropertyMap>
void thin_list_graph_with_priority(
  const PriorityValuePropertyMap& priority_value_property_map, const IsthmusProtectedPropertyMap& isthmus_protected_property_map,
  const ProtectedPropertyMap& protected_vertices,
  internal::Construction_feature_graph& graph, // [TODO] template, but check that clear_vertex and remove_vertex does not invalidate descriptors
  NextToHolePropertyMap& next_to_hole_property_map)
{
  typedef internal::Construction_feature_graph Graph;
  using Vertex = typename Graph::vertex_descriptor;
  using Sharpness_type = typename boost::property_traits<PriorityValuePropertyMap>::value_type;

  // sort vertices by priority values, and create lists for elements with the same values.
  std::map<Sharpness_type, std::list<Vertex>> sorted_elements;
  for (Vertex vertex : boost::make_iterator_range(vertices(graph)))
    sorted_elements[get(priority_value_property_map, vertex)].push_back(vertex);
  typename std::map<Sharpness_type, std::list<Vertex>>::iterator lowest_priority_to_check_iterator = sorted_elements.begin();

  // cache the removable property to not compute it multiple times
  std::vector<bool> is_element_computed_cache(num_vertices(graph), false);
  auto is_element_computed_cache_map = boost::make_iterator_property_map(is_element_computed_cache.begin(),
                                            get(boost::vertex_index, graph));

  bool deleted_simple = true;
  while (deleted_simple)
  {
    deleted_simple = internal::thin_remove_lowest_value_vertex(
      /*In    */ priority_value_property_map, isthmus_protected_property_map, protected_vertices,
      /*In|Out*/ graph, next_to_hole_property_map, sorted_elements, lowest_priority_to_check_iterator, is_element_computed_cache_map);
  }
}

} /* namespace Feature_graph */

} /* namespace CGAL */


#endif // CGAL_FEATURE_GRAPH_INTERNAL_GRAPH_OPERATION_H