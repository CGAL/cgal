#ifndef CGAL_SURFACE_MESH_SEGMENTATION_ALPHA_EXPANSION_GRAPH_CUT_H
// Copyright (c) 2014  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Ilker O. Yaz

#define CGAL_SURFACE_MESH_SEGMENTATION_ALPHA_EXPANSION_GRAPH_CUT_H

#include <CGAL/license/Surface_mesh_segmentation.h>


/// @cond CGAL_DOCUMENT_INTERNAL

/**
 * @file Alpha_expansion_graph_cut.h
 * @brief This file contains 3 graph-cut algorithms, which can be used as a template parameter for CGAL::internal::Surface_mesh_segmentation.
 *
 * Main differences between implementations are underlying max-flow algorithm and graph type (i.e. results are the same, performance differs).
 *
 * By default, we use MAXFLOW and the class Alpha_expansion_graph_cut_boykov_kolmogorov.
 * For deactivating MAXFLOW software and using boost implementation instead, define CGAL_DO_NOT_USE_BOYKOV_KOLMOGOROV_MAXFLOW_SOFTWARE.
 * It deactivates Alpha_expansion_graph_cut_boykov_kolmogorov, activate boost versions
 * and makes CGAL::internal::Surface_mesh_segmentation using Alpha_expansion_graph_cut_boost
 * as default implementation for the graph-cut.
 *
 * Also algorithms can be used by their-own for applying alpha-expansion graph-cut on any graph.
 *
 */
#include <CGAL/assertions.h>
#ifdef CGAL_SEGMENTATION_BENCH_GRAPHCUT
#include <CGAL/Timer.h>
#endif
#include <CGAL/IO/trace.h>

#include <CGAL/boost/graph/named_function_params.h>

#include <boost/version.hpp>
#ifdef CGAL_DO_NOT_USE_BOYKOV_KOLMOGOROV_MAXFLOW_SOFTWARE
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/compressed_sparse_row_graph.hpp>
#if BOOST_VERSION >= 104400 // at this version kolmogorov_max_flow become depricated.
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#else
#include <boost/graph/kolmogorov_max_flow.hpp>
#endif
#else
namespace MaxFlow
{
#include <CGAL/internal/auxiliary/graph.h>
}
#endif

#include <vector>




namespace CGAL
{
namespace internal
{

struct Alpha_expansion_old_API_wrapper_graph
{
  typedef std::size_t vertex_descriptor;
  typedef std::size_t edge_descriptor;
  typedef boost::directed_tag directed_category;
  typedef boost::disallow_parallel_edge_tag edge_parallel_category;
  typedef boost::edge_list_graph_tag traversal_category;
    
  typedef boost::counting_iterator<std::size_t> counting_iterator;
  typedef CGAL::Iterator_range<counting_iterator> counting_range;

  typedef CGAL::Identity_property_map<std::size_t> Vertex_index_map;
  typedef CGAL::Pointer_property_map<std::size_t>::type Vertex_label_map;

  struct Vertex_label_probability_map
  {
    typedef std::size_t key_type;
    typedef std::vector<double> value_type;
    typedef value_type reference;
    typedef boost::readable_property_map_tag category;
      
    const std::vector<std::vector<double> >* probability_matrix;

    Vertex_label_probability_map (const std::vector<std::vector<double> >* probability_matrix)
      : probability_matrix (probability_matrix)
    { }

    friend reference get (const Vertex_label_probability_map& pmap, key_type idx)
    {
      std::vector<double> out;
      out.reserve (pmap.probability_matrix->size());
      for (std::size_t i = 0; i < pmap.probability_matrix->size(); ++ i)
        out.push_back ((*pmap.probability_matrix)[i][idx]);
      return out;
    }

  };

  typedef CGAL::Pointer_property_map<double>::const_type Edge_weight_map;
    
  const std::vector<std::pair<std::size_t, std::size_t> >& edges;
  const std::vector<double>& edge_weights;
  const std::vector<std::vector<double> >& probability_matrix;
  std::vector<std::size_t>& labels;

  Alpha_expansion_old_API_wrapper_graph (const std::vector<std::pair<std::size_t, std::size_t> >& edges,
                                         const std::vector<double>& edge_weights,
                                         const std::vector<std::vector<double> >& probability_matrix,
                                         std::vector<std::size_t>& labels)
    : edges (edges), edge_weights (edge_weights), probability_matrix (probability_matrix), labels (labels)
  { }

  friend counting_range vertices (const Alpha_expansion_old_API_wrapper_graph& graph)
  {
    return CGAL::make_range (boost::counting_iterator<std::size_t>(0),
                             boost::counting_iterator<std::size_t>(graph.labels.size()));
  }

  friend std::size_t num_vertices (const Alpha_expansion_old_API_wrapper_graph& graph) { return graph.labels.size(); }

  friend counting_range edges (const Alpha_expansion_old_API_wrapper_graph& graph)
  {
    return CGAL::make_range (boost::counting_iterator<std::size_t>(0),
                             boost::counting_iterator<std::size_t>(graph.edges.size()));
  }

  friend vertex_descriptor source (edge_descriptor ed, const Alpha_expansion_old_API_wrapper_graph& graph)
  { return graph.edges[ed].first; }
  friend vertex_descriptor target (edge_descriptor ed, const Alpha_expansion_old_API_wrapper_graph& graph)
  { return graph.edges[ed].second; }

  Vertex_index_map vertex_index_map() const { return Vertex_index_map(); }
  Vertex_label_map vertex_label_map() { return CGAL::make_property_map(labels); }
  Vertex_label_probability_map vertex_label_probability_map() const
  { return Vertex_label_probability_map(&probability_matrix); }
  Edge_weight_map edge_weight_map() const { return CGAL::make_property_map(edge_weights); }
};

////////////////////////////////////////////////////////////////////////////////////////
//   Comments about performance:
//
// 1) With BGL:
//     * Using adjacency_list:
//     ** Without pre-allocating vertex-list
//       | OutEdgeList | VertexList | Performance |
//       |    listS    |   listS    |     25.2    |
//       |    vecS     |   listS    |     22.7    |
//       |    listS    |   vecS     |     30.7    |
//       |    vecS     |   vecS     |     26.1    |
//
//     ** With pre-allocating vertex-list with max-node size
//        (Note: exact number of vertices are not certain at the beginning)
//       | OutEdgeList | VertexList | Performance |
//       |    listS    |   vecS     |     25.2    |
//       |    vecS     |   vecS     |     23.4    |
//
//     * Didn't try adjacency_matrix since our graph is sparse
//     ( Also one can check BGL book, performance section )
//
//    Decision:
//     * Alpha_expansion_graph_cut_boost: use adjacency_list<vecS, listS> without
//       pre-allocating vertex-list.
//
// 2) With Boykov-Kolmogorov MAXFLOW software:
//   (http://pub.ist.ac.at/~vnk/software/maxflow-v2.21.src.tar.gz)
//                                  | Performance |
//                                  |     3.1     |
//     * Alpha_expansion_graph_cut_boykov_kolmogorov provides an implementation.
//       MAXFLOW does not provide any option for pre-allocation (It is possible with v_3.02 though).
//
// Typical Benchmark result provided by Ilker
//                                 | construction of vertices  |  construction of edges    | graph cut  | Total
//   -----------------------------------------------------------------------------------------------------------
//   boost with an adjacency list  |         1.53              | 1.51                      |  3.00      | 6.04
//   boost with CSR                | 0.11 (gather in a vector) | 0.15 (gather in a vector) |  2.67      | 2.93
//   MaxFlow                       |       0.042               | 0.076                     |  1.043     | 1.161
//
// The main issue for now with CSR is the construction of the opposite edge map that is too costly,
// since it is done by exploring all edges to find opposite
////////////////////////////////////////////////////////////////////////////////////////

#ifdef CGAL_DO_NOT_USE_BOYKOV_KOLMOGOROV_MAXFLOW_SOFTWARE
/**
 * @brief Implements alpha-expansion graph cut algorithm.
 *
 * For representing graph, it uses adjacency_list with OutEdgeList = vecS, VertexList = listS.
 * Also no pre-allocation is made for vertex-list.
 */
class Alpha_expansion_graph_cut_boost
{
private:
  typedef boost::adjacency_list_traits<boost::vecS, boost::listS, boost::directedS>
  Adjacency_list_traits;

  typedef boost::adjacency_list<boost::vecS, boost::listS, boost::directedS,
          // 4 vertex properties
          boost::property<boost::vertex_index_t, std::size_t,
          boost::property<boost::vertex_color_t, boost::default_color_type,
          boost::property<boost::vertex_distance_t, double,
          boost::property<boost::vertex_predecessor_t, Adjacency_list_traits::edge_descriptor >
          > > >,
          // 3 edge properties
          boost::property<boost::edge_capacity_t, double,
          boost::property<boost::edge_residual_capacity_t, double,
          boost::property<boost::edge_reverse_t, Adjacency_list_traits::edge_descriptor> >
          > > Graph;

  typedef boost::graph_traits<Graph> Traits;
  typedef boost::color_traits<boost::default_color_type> ColorTraits;

  typedef Traits::vertex_descriptor Vertex_descriptor;
  typedef Traits::vertex_iterator   Vertex_iterator;
  typedef Traits::edge_descriptor   Edge_descriptor;
  typedef Traits::edge_iterator     Edge_iterator;

  /**
   * Adds two directional edges between @a v1 and @a v2
   * @param v1 first vertex
   * @param v2 second vertex
   * @param w1 weight for edge from v1 to v2 (v1->v2)
   * @param w2 weight for edge from v2 to v1 (v2->v1)
   * @param graph to be added
   * @return pair of added edges, first: v1->v2 and second: v2->v1
   */
  boost::tuple<Edge_descriptor, Edge_descriptor>
  add_edge_and_reverse(Vertex_descriptor& v1, Vertex_descriptor& v2, double w1,
                       double w2, Graph& graph) const {
    Edge_descriptor v1_v2, v2_v1;
    bool v1_v2_added, v2_v1_added;

    boost::tie(v1_v2, v1_v2_added) = boost::add_edge(v1, v2, graph);
    boost::tie(v2_v1, v2_v1_added) = boost::add_edge(v2, v1, graph);

    CGAL_assertion(v1_v2_added && v2_v1_added);
    //put edge capacities
    boost::put(boost::edge_reverse, graph, v1_v2, v2_v1);
    boost::put(boost::edge_reverse, graph, v2_v1, v1_v2);

    //map reverse edges
    boost::put(boost::edge_capacity, graph, v1_v2, w1);
    boost::put(boost::edge_capacity, graph, v2_v1, w2);

    return boost::make_tuple(v1_v2, v2_v1);
  }



public:
  /**
   * Applies alpha-expansion graph-cut for energy minimization.
   * @param edges contains incident vertex-id pairs for each edge (vertex-ids should be between [0, number of vertices -1])
   * @param edge_weights contains weights for each edge in @a edges (correspondence according to order)
   * @param probability_matrix contains responsibility of the center on the vertex probability[center][vertex]
   * @param[in, out] labels as input it contains initial labeling of vertices (i.e. a center-id between [0, number of centers -1]),
   * and as output it returns final labeling of vertices (i.e. assigned cluster-id to each facet)
   * @return result of energy function
   */
  double operator()(const std::vector<std::pair<std::size_t, std::size_t> >&
                    edges,
                    const std::vector<double>& edge_weights,
                    const std::vector<std::vector<double> >& probability_matrix,
                    std::vector<std::size_t>& labels) const
  {
    Alpha_expansion_old_API_wrapper_graph graph (edges, edge_weights, probability_matrix, labels);

    return (*this)(graph,
                   graph.edge_weight_map(),
                   graph.vertex_index_map(),
                   graph.vertex_label_map(),
                   graph.vertex_label_probability_map());
  }

  template <typename InputGraph,
            typename Edge_weight_map,
            typename Vertex_index_map,
            typename Vertex_label_map,
            typename Vertex_label_probability_map>
  double operator()(const InputGraph& input_graph,
                    Edge_weight_map edge_weight_map,
                    Vertex_index_map vertex_index_map,
                    Vertex_label_map vertex_label_map,
                    Vertex_label_probability_map vertex_label_probability_map) const
  {
    typedef boost::graph_traits<InputGraph> GT;
    typedef typename GT::edge_descriptor input_edge_descriptor;
    typedef typename GT::vertex_descriptor input_vertex_descriptor;

    // TODO: check this hardcoded parameter
    const double tolerance = 1e-10;

    double min_cut = (std::numeric_limits<double>::max)();

    #ifdef CGAL_SEGMENTATION_BENCH_GRAPHCUT
    double vertex_creation_time, edge_creation_time, cut_time;
    vertex_creation_time = edge_creation_time = cut_time = 0.0;
    #endif

    std::vector<Vertex_descriptor> inserted_vertices;
    inserted_vertices.resize(num_vertices (input_graph));

    std::size_t number_of_labels = get(vertex_label_probability_map, *std::begin(vertices(input_graph))).size();

    Graph graph;

    bool success;
    do {
      success = false;

      for (std::size_t alpha = 0; alpha < number_of_labels; ++ alpha)
      {
        graph.clear();

        Vertex_descriptor cluster_source = boost::add_vertex(graph);
        Vertex_descriptor cluster_sink = boost::add_vertex(graph);

        #ifdef CGAL_SEGMENTATION_BENCH_GRAPHCUT
        Timer timer;
        timer.start();
        #endif

        // For E-Data
        // add every input vertex as a vertex to the graph, put edges to source & sink vertices
        for (input_vertex_descriptor vd : vertices(input_graph))
        {
          std::size_t vertex_i = get(vertex_index_map, vd);
          Vertex_descriptor new_vertex = boost::add_vertex(graph);
          inserted_vertices[vertex_i] = new_vertex;
          double source_weight = get(vertex_label_probability_map, vd)[alpha];
          // since it is expansion move, current alpha labeled vertices will be assigned to alpha again,
          // making sink_weight 'infinity' guarantee this.
          double sink_weight = (get(vertex_label_map, vd) == alpha ?
                                (std::numeric_limits<double>::max)()
                                : get(vertex_label_probability_map, vd)[get(vertex_label_map, vd)]);

          add_edge_and_reverse(cluster_source, new_vertex, source_weight, 0.0, graph);
          add_edge_and_reverse(new_vertex, cluster_sink, sink_weight, 0.0, graph);
        }
        #ifdef CGAL_SEGMENTATION_BENCH_GRAPHCUT
        vertex_creation_time += timer.time();
        timer.reset();
        #endif

        // For E-Smooth
        // add edge between every vertex,
        for (input_edge_descriptor ed : edges(input_graph))
        {
          input_vertex_descriptor vd1 = source(ed, input_graph);
          input_vertex_descriptor vd2 = target(ed, input_graph);
          std::size_t idx1 = get (vertex_index_map, vd1);
          std::size_t idx2 = get (vertex_index_map, vd2);
          
          double weight = get (edge_weight_map, ed);
          
          Vertex_descriptor v1 = inserted_vertices[idx1],
                            v2 = inserted_vertices[idx2];
          
          std::size_t label_1 = get (vertex_label_map, vd1);
          std::size_t label_2 = get (vertex_label_map, vd2);
          if(label_1 == label_2) {
            if(label_1 != alpha) {
              add_edge_and_reverse(v1, v2, weight, weight, graph);
            }
          } else {
            Vertex_descriptor inbetween = boost::add_vertex(graph);

            double w1 = (label_1 == alpha) ? 0 : weight;
            double w2 = (label_2 == alpha) ? 0 : weight;
            add_edge_and_reverse(inbetween, v1, w1, w1, graph);
            add_edge_and_reverse(inbetween, v2, w2, w2, graph);
            add_edge_and_reverse(inbetween, cluster_sink, weight, 0.0, graph);
          }
        }
        #ifdef CGAL_SEGMENTATION_BENCH_GRAPHCUT
        edge_creation_time += timer.time();
        #endif

        // initialize vertex indices, it is necessary since we are using VertexList = listS
        Vertex_iterator v_begin, v_end;
        Traits::vertices_size_type index = 0;
        for(boost::tie(v_begin, v_end) = vertices(graph); v_begin != v_end; ++v_begin) {
          boost::put(boost::vertex_index, graph, *v_begin, index++);
        }

        #ifdef CGAL_SEGMENTATION_BENCH_GRAPHCUT
        timer.reset();
        #endif
#if BOOST_VERSION >= 104400
        double flow = boost::boykov_kolmogorov_max_flow(graph, cluster_source,
                      cluster_sink);
#else
        double flow = boost::kolmogorov_max_flow(graph, cluster_source, cluster_sink);
#endif
        #ifdef CGAL_SEGMENTATION_BENCH_GRAPHCUT
        cut_time += timer.time();
        #endif

        if(min_cut - flow <= flow * tolerance) {
          continue;
        }
        min_cut = flow;
        success = true;
        //update labeling
        for (input_vertex_descriptor vd : vertices (input_graph))
        {
          std::size_t vertex_i = get (vertex_index_map, vd);
          boost::default_color_type color = boost::get(boost::vertex_color, graph,
                                            inserted_vertices[vertex_i]);
          if(get (vertex_label_map, vd) != alpha
              && color == ColorTraits::white()) { //new comers (expansion occurs)
            put (vertex_label_map, vd, alpha);
          }
        }
      }
    } while(success);

    #ifdef CGAL_SEGMENTATION_BENCH_GRAPHCUT
    CGAL_TRACE_STREAM << "vertex creation time: " << vertex_creation_time <<
                      std::endl;
    CGAL_TRACE_STREAM << "edge creation time: " << edge_creation_time << std::endl;
    CGAL_TRACE_STREAM << "max flow algorithm time: " << cut_time << std::endl;
    #endif

    return min_cut;
  }
};

// another implementation using compressed_sparse_row_graph
// for now there is a performance problem while setting reverse edges
// if that can be solved, it is faster than Alpha_expansion_graph_cut_boost
class Alpha_expansion_graph_cut_boost_CSR
{
private:
  // CSR only accepts bundled props
  struct VertexP {
    boost::default_color_type vertex_color;
    double vertex_distance_t;
    // ? do not now there is another way to take it, I think since edge_descriptor does not rely on properties
    // this should be fine...
    boost::compressed_sparse_row_graph<boost::directedS>::edge_descriptor
    vertex_predecessor;
  };

  struct EdgeP {
    double edge_capacity;
    double edge_residual_capacity;
    boost::compressed_sparse_row_graph<boost::directedS>::edge_descriptor
    edge_reverse;
  };

  typedef boost::compressed_sparse_row_graph<boost::directedS,
          VertexP, EdgeP> Graph;

  typedef boost::graph_traits<Graph> Traits;
  typedef boost::color_traits<boost::default_color_type> ColorTraits;

  typedef Traits::vertex_descriptor Vertex_descriptor;
  typedef Traits::vertex_iterator   Vertex_iterator;
  typedef Traits::edge_descriptor   Edge_descriptor;
  typedef Traits::edge_iterator     Edge_iterator;

  void
  add_edge_and_reverse(std::size_t v1 , std::size_t v2, double w1, double w2,
                       std::vector<std::pair<std::size_t, std::size_t> >& edge_map,
                       std::vector<EdgeP>& edge_weights) const {
    edge_map.push_back(std::make_pair(v1, v2));
    EdgeP p1;
    p1.edge_capacity = w1;
    edge_weights.push_back(p1);

    edge_map.push_back(std::make_pair(v2, v1));
    EdgeP p2;
    p2.edge_capacity = w2;
    edge_weights.push_back(p2);
  }

public:
  /**
   * Applies alpha-expansion graph-cut for energy minimization.
   * @param edges contains incident vertex-id pairs for each edge (vertex-ids should be between [0, number of vertices -1])
   * @param edge_weights contains weights for each edge in @a edges (correspondence according to order)
   * @param probability_matrix contains responsibility of the center on the vertex probability[center][vertex]
   * @param[in, out] labels as input it contains initial labeling of vertices (i.e. a center-id between [0, number of centers -1]),
   * and as output it returns final labeling of vertices (i.e. assigned cluster-id to each facet)
   * @return result of energy function
   */
  double operator()(const std::vector<std::pair<std::size_t, std::size_t> >&
                    edges,
                    const std::vector<double>& edge_weights,
                    const std::vector<std::vector<double> >& probability_matrix,
                    std::vector<std::size_t>& labels) const {
    const double tolerance = 1e-10;

    double min_cut = (std::numeric_limits<double>::max)();

    #ifdef CGAL_SEGMENTATION_BENCH_GRAPHCUT
    double vertex_creation_time, edge_creation_time, graph_creation_time,
           reverse_mapping_time, cut_time;
    vertex_creation_time = edge_creation_time = graph_creation_time =
                             reverse_mapping_time = cut_time = 0.0;
    #endif

    Graph graph;

    bool success;
    do {
      success = false;
      std::size_t alpha = 0;

      for(std::vector<std::vector<double> >::const_iterator it =
            probability_matrix.begin();
          it != probability_matrix.end(); ++it, ++alpha) {
        std::vector<std::pair<std::size_t, std::size_t> > edge_map;
        std::vector<EdgeP>                edge_map_weights;
        edge_map.reserve(labels.size() *
                         8); // there is no way to know exact edge count, it is a heuristic value
        edge_map_weights.reserve(labels.size() * 8);

        std::size_t cluster_source = 0;
        std::size_t cluster_sink = 1;

        #ifdef CGAL_SEGMENTATION_BENCH_GRAPHCUT
        Timer timer;
        timer.start();
        #endif
        // For E-Data
        // add every facet as a vertex to the graph, put edges to source & sink vertices
        for(std::size_t vertex_i = 0; vertex_i < labels.size(); ++vertex_i) {
          double source_weight = probability_matrix[alpha][vertex_i];
          // since it is expansion move, current alpha labeled vertices will be assigned to alpha again,
          // making sink_weight 'infinity' guarantee this.
          double sink_weight = (labels[vertex_i] == alpha) ?
                               (std::numeric_limits<double>::max)()
                               : probability_matrix[labels[vertex_i]][vertex_i];

          add_edge_and_reverse(cluster_source, vertex_i + 2, source_weight, 0.0, edge_map,
                               edge_map_weights);
          add_edge_and_reverse(vertex_i + 2, cluster_sink, sink_weight, 0.0, edge_map,
                               edge_map_weights);
        }
        #ifdef CGAL_SEGMENTATION_BENCH_GRAPHCUT
        vertex_creation_time += timer.time();
        timer.reset();
        #endif
        // For E-Smooth
        // add edge between every vertex,
        std::size_t num_vert = labels.size() + 2;
        std::vector<double>::const_iterator weight_it = edge_weights.begin();
        for(std::vector<std::pair<std::size_t, std::size_t> >::const_iterator edge_it =
              edges.begin(); edge_it != edges.end();
            ++edge_it, ++weight_it) {
          std::size_t v1 = edge_it->first + 2, v2 = edge_it->second + 2;
          std::size_t label_1 = labels[edge_it->first], label_2 = labels[edge_it->second];
          if(label_1 == label_2) {
            if(label_1 != alpha) {
              add_edge_and_reverse(v1, v2, *weight_it, *weight_it, edge_map,
                                   edge_map_weights);
            }
          } else {
            std::size_t inbetween = num_vert++;

            double w1 = (label_1 == alpha) ? 0 : *weight_it;
            double w2 = (label_2 == alpha) ? 0 : *weight_it;
            add_edge_and_reverse(inbetween, v1, w1, w1, edge_map, edge_map_weights);
            add_edge_and_reverse(inbetween, v2, w2, w2, edge_map, edge_map_weights);
            add_edge_and_reverse(inbetween, cluster_sink, *weight_it, 0.0, edge_map,
                                 edge_map_weights);
          }
        }
        #ifdef CGAL_SEGMENTATION_BENCH_GRAPHCUT
        edge_creation_time += timer.time();
        timer.reset();
        #endif
#if BOOST_VERSION >= 104000
        Graph graph(boost::edges_are_unsorted, edge_map.begin(), edge_map.end(),
                    edge_map_weights.begin(), num_vert);
#else
        Graph graph(edge_map.begin(), edge_map.end(),
                    edge_map_weights.begin(), num_vert);
#endif
        #ifdef CGAL_SEGMENTATION_BENCH_GRAPHCUT
        graph_creation_time += timer.time();
        timer.reset();
        #endif

        // PERFORMANCE PROBLEM
        // need to set reverse edge map, I guess there is no way to do that before creating the graph
        // since we do not have edge_descs
        // however from our edge_map, we know that each (2i, 2i + 1) is reverse pairs, how to facilitate that ?
        // will look it back
        Graph::edge_iterator ei, ee;
        for(boost::tie(ei, ee) = boost::edges(graph); ei != ee; ++ei) {
          Graph::vertex_descriptor v1 = boost::source(*ei, graph);
          Graph::vertex_descriptor v2 = boost::target(*ei, graph);
          std::pair<Graph::edge_descriptor, bool> opp_edge = boost::edge(v2, v1, graph);

          CGAL_assertion(opp_edge.second);
          graph[opp_edge.first].edge_reverse =
            *ei; // and edge_reverse of *ei will be (or already have been) set by the opp_edge
        }
        #ifdef CGAL_SEGMENTATION_BENCH_GRAPHCUT
        reverse_mapping_time += timer.time();
        timer.reset();
        #endif

#if BOOST_VERSION >= 104400
        // since properties are bundled, defaults does not work need to specify them
        double flow = boost::boykov_kolmogorov_max_flow(graph,
                      boost::get(&EdgeP::edge_capacity, graph),
                      boost::get(&EdgeP::edge_residual_capacity, graph),
                      boost::get(&EdgeP::edge_reverse, graph),
                      boost::get(&VertexP::vertex_predecessor, graph),
                      boost::get(&VertexP::vertex_color, graph),
                      boost::get(&VertexP::vertex_distance_t, graph),
                      boost::get(boost::vertex_index,
                                 graph), // this is not bundled, get it from graph (CRS provides one)
                      cluster_source,
                      cluster_sink
                                                       );
#else
        double flow = boost::kolmogorov_max_flow(graph,
                      boost::get(&EdgeP::edge_capacity, graph),
                      boost::get(&EdgeP::edge_residual_capacity, graph),
                      boost::get(&EdgeP::edge_reverse, graph),
                      boost::get(&VertexP::vertex_predecessor, graph),
                      boost::get(&VertexP::vertex_color, graph),
                      boost::get(&VertexP::vertex_distance_t, graph),
                      boost::get(boost::vertex_index,
                                 graph), // this is not bundled, get it from graph
                      cluster_source,
                      cluster_sink
                                                );
#endif
        #ifdef CGAL_SEGMENTATION_BENCH_GRAPHCUT
        cut_time += timer.time();
        #endif

        if(min_cut - flow <= flow * tolerance) {
          continue;
        }
        min_cut = flow;
        success = true;
        //update labeling
        for(std::size_t vertex_i = 0; vertex_i < labels.size(); ++vertex_i) {
          boost::default_color_type color =  graph[vertex_i + 2].vertex_color;
          if(labels[vertex_i] != alpha
              && color == ColorTraits::white()) { //new comers (expansion occurs)
            labels[vertex_i] = alpha;
          }
        }
      }
    } while(success);

    #ifdef CGAL_SEGMENTATION_BENCH_GRAPHCUT
    CGAL_TRACE_STREAM << "vertex creation time: " << vertex_creation_time <<
                      std::endl;
    CGAL_TRACE_STREAM << "edge creation time: " << edge_creation_time << std::endl;
    CGAL_TRACE_STREAM << "graph creation time: " << graph_creation_time <<
                      std::endl;
    CGAL_TRACE_STREAM << "reverse mapping time: " << reverse_mapping_time <<
                      std::endl;
    CGAL_TRACE_STREAM << "max flow algorithm time: " << cut_time << std::endl;
    #endif
    return min_cut;
  }
};
#endif

#ifndef CGAL_DO_NOT_USE_BOYKOV_KOLMOGOROV_MAXFLOW_SOFTWARE
/**
 * @brief Implements alpha-expansion graph cut algorithm.
 *
 * For underlying max-flow algorithm, it uses the MAXFLOW software implemented by Boykov & Kolmogorov.
 *  Also no pre-allocation is made.
 */
class Alpha_expansion_graph_cut_boykov_kolmogorov
{
public:
  /**
   * Applies alpha-expansion graph-cut for energy minimization.
   * @param edges contains incident vertex-id pairs for each edge (vertex-ids should be between [0, number of vertices -1])
   * @param edge_weights contains weights for each edge in @a edges (correspondence according to order)
   * @param probability_matrix contains responsibility of the center on the vertex probability[center][vertex]
   * @param[in, out] labels as input it contains initial labeling of vertices (i.e. a center-id between [0, number of centers -1]),
   * and as output it returns final labeling of vertices
   * @return result of energy function
   */
  double operator()(const std::vector<std::pair<std::size_t, std::size_t> >&
                    edges,
                    const std::vector<double>& edge_weights,
                    const std::vector<std::vector<double> >& probability_matrix,
                    std::vector<std::size_t>& labels) const {
    const double tolerance = 1e-10;

    double min_cut = (std::numeric_limits<double>::max)();

    #ifdef CGAL_SEGMENTATION_BENCH_GRAPHCUT
    double vertex_creation_time, edge_creation_time, cut_time;
    vertex_creation_time = edge_creation_time = cut_time = 0.0;
    #endif

    std::vector<MaxFlow::Graph::node_id> inserted_vertices;
    inserted_vertices.resize(labels.size());
    bool success;
    do {
      success = false;
      std::size_t alpha = 0;
      for(std::vector<std::vector<double> >::const_iterator it =
            probability_matrix.begin();
          it != probability_matrix.end(); ++it, ++alpha) {
        MaxFlow::Graph graph;
        #ifdef CGAL_SEGMENTATION_BENCH_GRAPHCUT
        Timer timer;
        timer.start();
        #endif
        // For E-Data
        // add every facet as a vertex to graph, put edges to source-sink vertices
        for(std::size_t vertex_i = 0; vertex_i <  labels.size(); ++vertex_i) {
          MaxFlow::Graph::node_id new_vertex = graph.add_node();
          inserted_vertices[vertex_i] = new_vertex;

          double source_weight = probability_matrix[alpha][vertex_i];
          // since it is expansion move, current alpha labeled vertices will be assigned to alpha again,
          // making sink_weight 'infinity' guarantee this.
          double sink_weight = (labels[vertex_i] == alpha) ?
                               (std::numeric_limits<double>::max)()
                               : probability_matrix[labels[vertex_i]][vertex_i];
          graph.add_tweights(new_vertex, source_weight, sink_weight);
        }
        #ifdef CGAL_SEGMENTATION_BENCH_GRAPHCUT
        vertex_creation_time += timer.time();
        timer.reset();
        #endif
        // For E-Smooth
        // add edge between every vertex,
        std::vector<double>::const_iterator weight_it = edge_weights.begin();
        for(std::vector<std::pair<std::size_t, std::size_t> >::const_iterator edge_it =
              edges.begin(); edge_it != edges.end();
            ++edge_it, ++weight_it) {
          MaxFlow::Graph::node_id v1 = inserted_vertices[edge_it->first];
          MaxFlow::Graph::node_id v2 = inserted_vertices[edge_it->second];
          std::size_t label_1 = labels[edge_it->first], label_2 = labels[edge_it->second];
          if(label_1 == label_2) {
            if(label_1 != alpha) {
              graph.add_edge(v1, v2, *weight_it, *weight_it);
            }
          } else {
            MaxFlow::Graph::node_id inbetween = graph.add_node();

            double w1 = (label_1 == alpha) ? 0 : *weight_it;
            double w2 = (label_2 == alpha) ? 0 : *weight_it;
            graph.add_edge(inbetween, v1, w1, w1);
            graph.add_edge(inbetween, v2, w2, w2);

            graph.add_tweights(inbetween, 0.0, *weight_it);
          }
        }
        #ifdef CGAL_SEGMENTATION_BENCH_GRAPHCUT
        edge_creation_time += timer.time();
        timer.reset();
        #endif

        double flow = graph.maxflow();
        #ifdef CGAL_SEGMENTATION_BENCH_GRAPHCUT
        cut_time += timer.time();
        #endif

        if(min_cut - flow <= flow * tolerance) {
          continue;
        }

        min_cut = flow;
        success = true;
        //update labeling
        for(std::size_t vertex_i = 0; vertex_i < labels.size(); ++vertex_i) {
          if(labels[vertex_i] != alpha
              && graph.what_segment(inserted_vertices[vertex_i]) == MaxFlow::Graph::SINK) {
            labels[vertex_i] = alpha;
          }
        }
      }
    } while(success);

    #ifdef CGAL_SEGMENTATION_BENCH_GRAPHCUT
    CGAL_TRACE_STREAM << "vertex creation time: " << vertex_creation_time <<
                      std::endl;
    CGAL_TRACE_STREAM << "edge creation time: " << edge_creation_time << std::endl;
    CGAL_TRACE_STREAM << "max flow algorithm time: " << cut_time << std::endl;
    #endif
    return min_cut;
  }
};
#endif //CGAL_DO_NOT_USE_BOYKOV_KOLMOGOROV_MAXFLOW_SOFTWARE
}//namespace internal
/// @endcond

struct Alpha_expansion_boost_adjacency_list { };
struct Alpha_expansion_boost_compressed_sparse_row { };
struct Alpha_expansion_MaxFlow { };
  

template <typename InputGraph,
          typename Edge_weight_map,
          typename Vertex_index_map,
          typename Vertex_label_map,
          typename Vertex_label_probability_map,
          typename AlphaExpansionImplementation>
double alpha_expansion_graph_cut (const InputGraph& input_graph,
                                  Edge_weight_map edge_weight_map,
                                  Vertex_index_map vertex_index_map,
                                  Vertex_label_map vertex_label_map,
                                  Vertex_label_probability_map vertex_label_probability_map,
                                  const AlphaExpansionImplementation& = AlphaExpansionImplementation())
{
  if (std::is_same<AlphaExpansionImplementation,  Alpha_expansion_boost_adjacency_list>::value)
    return internal::Alpha_expansion_graph_cut_boost()(input_graph,
                                                       edge_weight_map,
                                                       vertex_index_map,
                                                       vertex_label_map,
                                                       vertex_label_probability_map);
}

template <typename InputGraph,
          typename Edge_weight_map,
          typename Vertex_index_map,
          typename Vertex_label_map,
          typename Vertex_label_probability_map>
double alpha_expansion_graph_cut (const InputGraph& input_graph,
                                  Edge_weight_map edge_weight_map,
                                  Vertex_index_map vertex_index_map,
                                  Vertex_label_map vertex_label_map,
                                  Vertex_label_probability_map vertex_label_probability_map)
{
  return alpha_expansion_graph_cut<InputGraph,
                                   Edge_weight_map,
                                   Vertex_index_map,
                                   Vertex_label_map,
                                   Vertex_label_probability_map,
                                   Alpha_expansion_boost_adjacency_list>
    (input_graph, edge_weight_map, vertex_index_map, vertex_label_map, vertex_label_probability_map);
}

}//namespace CGAL
#endif //CGAL_SURFACE_MESH_SEGMENTATION_ALPHA_EXPANSION_GRAPH_CUT_H
