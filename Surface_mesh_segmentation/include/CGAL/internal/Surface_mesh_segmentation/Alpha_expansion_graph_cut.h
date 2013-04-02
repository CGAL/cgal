#ifndef CGAL_SURFACE_MESH_SEGMENTATION_ALPHA_EXPANSION_GRAPH_CUT_H
#define CGAL_SURFACE_MESH_SEGMENTATION_ALPHA_EXPANSION_GRAPH_CUT_H

/// @cond CGAL_DOCUMENT_INTERNAL

/**
 * @file Alpha_expansion_graph_cut.h
 * @brief This file contains 3 graph-cut algorithms, which can be used as a template parameter for CGAL::internal::Surface_mesh_segmentation.
 *
 * Main differences between implementations are underlying max-flow algorithm and graph type (i.e. results are the same, performance differs).
 *
 * For activating MAXFLOW software, define CGAL_USE_BOYKOV_KOLMOGOROV_MAXFLOW_SOFTWARE. It activates Alpha_expansion_graph_cut_boykov_kolmogorov,
 * and makes CGAL::internal::Surface_mesh_segmentation choose this implementation for graph-cut.
 *
 * Also algorithms can be used by their-own for applying alpha-expansion graph-cut on any graph.
 */
#include <CGAL/assertions.h>

#include <boost/version.hpp>
#include <boost/graph/adjacency_list.hpp>
#if BOOST_VERSION >= 104400 // at this version kolmogorov_max_flow become depricated.
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#else
#include <boost/graph/kolmogorov_max_flow.hpp>
#endif

#include <vector>

//#define CGAL_USE_BOYKOV_KOLMOGOROV_MAXFLOW_SOFTWARE

#ifdef CGAL_USE_BOYKOV_KOLMOGOROV_MAXFLOW_SOFTWARE
#include <CGAL/internal/auxiliary/graph.h>
#endif

namespace CGAL
{
namespace internal
{

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
////////////////////////////////////////////////////////////////////////////////////////

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
  double operator()(const std::vector<std::pair<int, int> >& edges,
                    const std::vector<double>& edge_weights,
                    const std::vector<std::vector<double> >& probability_matrix,
                    std::vector<int>& labels) const {
    double min_cut = (std::numeric_limits<double>::max)();
    bool success;
    do {
      success = false;
      int alpha = 0;
      for(std::vector<std::vector<double> >::const_iterator it =
            probability_matrix.begin();
          it != probability_matrix.end(); ++it, ++alpha) {
        Graph graph;
        Vertex_descriptor cluster_source = boost::add_vertex(graph);
        Vertex_descriptor cluster_sink = boost::add_vertex(graph);
        std::vector<Vertex_descriptor> inserted_vertices;
        inserted_vertices.reserve(labels.size());

        // For E-Data
        // add every facet as a vertex to the graph, put edges to source & sink vertices
        for(std::size_t vertex_i = 0; vertex_i < labels.size(); ++vertex_i) {
          Vertex_descriptor new_vertex = boost::add_vertex(graph);
          inserted_vertices.push_back(new_vertex);
          double source_weight = probability_matrix[alpha][vertex_i];
          // since it is expansion move, current alpha labeled vertices will be assigned to alpha again,
          // making sink_weight 'infinity' guarantee this.
          double sink_weight = (labels[vertex_i] == alpha) ?
                               (std::numeric_limits<double>::max)()
                               : probability_matrix[labels[vertex_i]][vertex_i];

          add_edge_and_reverse(cluster_source, new_vertex, source_weight, 0.0, graph);
          add_edge_and_reverse(new_vertex, cluster_sink, sink_weight, 0.0, graph);
        }

        // For E-Smooth
        // add edge between every vertex,
        std::vector<double>::const_iterator weight_it = edge_weights.begin();
        for(std::vector<std::pair<int, int> >::const_iterator edge_it = edges.begin();
            edge_it != edges.end();
            ++edge_it, ++weight_it) {
          Vertex_descriptor v1 = inserted_vertices[edge_it->first],
                            v2 = inserted_vertices[edge_it->second];
          int label_1 = labels[edge_it->first], label_2 = labels[edge_it->second];
          if(label_1 == label_2) {
            if(label_1 != alpha) {
              add_edge_and_reverse(v1, v2, *weight_it, *weight_it, graph);
            }
          } else {
            Vertex_descriptor inbetween = boost::add_vertex(graph);

            double w1 = (label_1 == alpha) ? 0 : *weight_it;
            double w2 = (label_2 == alpha) ? 0 : *weight_it;
            add_edge_and_reverse(inbetween, v1, w1, w1, graph);
            add_edge_and_reverse(inbetween, v2, w2, w2, graph);
            add_edge_and_reverse(inbetween, cluster_sink, *weight_it, 0.0, graph);
          }
        }
        // initialize vertex indices, it is neccessary since we are using VertexList = listS
        Vertex_iterator v_begin, v_end;
        Traits::vertices_size_type index = 0;
        for(boost::tie(v_begin, v_end) = vertices(graph); v_begin != v_end; ++v_begin) {
          boost::put(boost::vertex_index, graph, *v_begin, index++);
        }

#if BOOST_VERSION >= 104400
        double flow = boost::boykov_kolmogorov_max_flow(graph, cluster_source,
                      cluster_sink);
#else
        double flow = boost::kolmogorov_max_flow(graph, cluster_source, cluster_sink);
#endif

        if(min_cut - flow < flow * 1e-10) {
          continue;
        }
        min_cut = flow;
        success = true;
        //update labeling
        for(std::size_t vertex_i = 0; vertex_i < inserted_vertices.size(); ++vertex_i) {
          boost::default_color_type color = boost::get(boost::vertex_color, graph,
                                            inserted_vertices[vertex_i]);
          if(labels[vertex_i] != alpha
              && color == ColorTraits::white()) { //new comers (expansion occurs)
            labels[vertex_i] = alpha;
          }
        }

      }
    } while(success);
    return min_cut;
  }
};

#ifdef CGAL_USE_BOYKOV_KOLMOGOROV_MAXFLOW_SOFTWARE
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
  double operator()(const std::vector<std::pair<int, int> >& edges,
                    const std::vector<double>& edge_weights,
                    const std::vector<std::vector<double> >& probability_matrix,
                    std::vector<int>& labels) const {
    double min_cut = (std::numeric_limits<double>::max)();
    bool success;
    do {
      success = false;
      int alpha = 0;
      for(std::vector<std::vector<double> >::const_iterator it =
            probability_matrix.begin();
          it != probability_matrix.end(); ++it, ++alpha) {
        Graph graph;
        std::vector<Graph::node_id> inserted_vertices;
        inserted_vertices.reserve(labels.size());
        // For E-Data
        // add every facet as a vertex to graph, put edges to source-sink vertices
        for(std::size_t vertex_i = 0; vertex_i <  probability_matrix[0].size();
            ++vertex_i) {
          Graph::node_id new_vertex = graph.add_node();
          inserted_vertices.push_back(new_vertex);

          double source_weight = probability_matrix[alpha][vertex_i];
          // since it is expansion move, current alpha labeled vertices will be assigned to alpha again,
          // making sink_weight 'infinity' guarantee this.
          double sink_weight = (labels[vertex_i] == alpha) ?
                               (std::numeric_limits<double>::max)()
                               : probability_matrix[labels[vertex_i]][vertex_i];
          graph.add_tweights(new_vertex, source_weight, sink_weight);
        }
        // For E-Smooth
        // add edge between every vertex,
        std::vector<double>::const_iterator weight_it = edge_weights.begin();
        for(std::vector<std::pair<int, int> >::const_iterator edge_it = edges.begin();
            edge_it != edges.end();
            ++edge_it, ++weight_it) {
          Graph::node_id v1 = inserted_vertices[edge_it->first];
          Graph::node_id v2 = inserted_vertices[edge_it->second];
          int label_1 = labels[edge_it->first], label_2 = labels[edge_it->second];
          if(label_1 == label_2) {
            if(label_1 != alpha) {
              graph.add_edge(v1, v2, *weight_it, *weight_it);
            }
          } else {
            Graph::node_id inbetween = graph.add_node();

            double w1 = (label_1 == alpha) ? 0 : *weight_it;
            double w2 = (label_2 == alpha) ? 0 : *weight_it;
            graph.add_edge(inbetween, v1, w1, w1);
            graph.add_edge(inbetween, v2, w2, w2);

            graph.add_tweights(inbetween, 0.0, *weight_it);
          }
        }
        double flow = graph.maxflow();
        if(min_cut - flow < flow * 1e-10) {
          continue;
        }

        min_cut = flow;
        success = true;
        //update labeling
        for(std::size_t vertex_i = 0; vertex_i < labels.size(); ++vertex_i) {
          if(labels[vertex_i] != alpha
              && graph.what_segment(inserted_vertices[vertex_i]) == Graph::SINK) {
            labels[vertex_i] = alpha;
          }
        }
      }
    } while(success);
    return min_cut;
  }
};
#endif //CGAL_USE_BOYKOV_KOLMOGOROV_MAXFLOW_SOFTWARE
}//namespace internal
/// @endcond
}//namespace CGAL
#endif //CGAL_SURFACE_MESH_SEGMENTATION_ALPHA_EXPANSION_GRAPH_CUT_H
