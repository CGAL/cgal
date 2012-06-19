#ifndef CGAL_ALPHA_EXPANSION_GRAPH_CUT_H
#define CGAL_ALPHA_EXPANSION_GRAPH_CUT_H

#include <CGAL/assertions.h>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#include <boost/math/distributions/normal.hpp>

#include <iostream>
#include <vector>

namespace CGAL
{
namespace internal
{

class Alpha_expansion_graph_cut
{
public:
  typedef boost::adjacency_list_traits<boost::vecS, boost::vecS, boost::directedS>
  Adjacency_list_traits;

  typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,
          // 4 vertex properties (nested)
          //boost::no_property
          boost::property<boost::vertex_index_t, int,
          boost::property<boost::vertex_color_t, boost::default_color_type,
          boost::property<boost::vertex_distance_t, double,
          boost::property<boost::vertex_predecessor_t, Adjacency_list_traits::edge_descriptor >
          > > >,
          // 3 edge properties (nested)
          boost::property<boost::edge_capacity_t, double,
          boost::property<boost::edge_residual_capacity_t, double,
          boost::property<boost::edge_reverse_t, Adjacency_list_traits::edge_descriptor> >
          > > Graph;

  typedef boost::graph_traits<Graph> Traits;
  typedef boost::color_traits<boost::default_color_type> ColorTraits;

  typedef Traits::vertex_descriptor Vertex_descriptor;
  typedef Traits::edge_descriptor Edge_descriptor;
  typedef Traits::vertex_iterator Vertex_iterator;
  typedef Traits::edge_iterator Edge_iterator;

  Alpha_expansion_graph_cut(const std::vector<std::pair<int, int> >& edges,
                            const std::vector<double>& edge_weights, std::vector<int>& labels,
                            const std::vector<std::vector<double> >& probability_matrix,
                            std::vector<int>& center_ids) {
    apply_alpha_expansion_2(edges, edge_weights, probability_matrix, labels);
    center_ids = labels;
  }

  boost::tuple<Edge_descriptor, Edge_descriptor>
  add_edge_and_reverse(Vertex_descriptor& v1, Vertex_descriptor& v2, double w1,
                       double w2, Graph& graph) {
    Edge_descriptor v1_v2, v2_v1;
    bool v1_v2_added, v2_v1_added;

    tie(v1_v2, v1_v2_added) = boost::add_edge(v1, v2, graph);
    tie(v2_v1, v2_v1_added) = boost::add_edge(v2, v1, graph);

    CGAL_assertion(v1_v2_added && v2_v1_added);
    //put edge capacities
    boost::put(boost::edge_reverse, graph, v1_v2, v2_v1);
    boost::put(boost::edge_reverse, graph, v2_v1, v1_v2);

    //map reverse edges
    boost::put(boost::edge_capacity, graph, v1_v2, w1);
    boost::put(boost::edge_capacity, graph, v2_v1, w2);

    return boost::make_tuple(v1_v2, v2_v1);
  }

  void apply_alpha_expansion(const std::vector<std::pair<int, int> >& edges,
                             const std::vector<double>& edge_weights,
                             const std::vector<std::vector<double> >& probability_matrix,
                             std::vector<int>& labels) {
    int number_of_clusters = probability_matrix.size();
    double min_cut = (std::numeric_limits<double>::max)();
    bool success;
    do {
      success = false;
      for(int alpha = 0; alpha < number_of_clusters; ++alpha) {
        Graph graph;
        Vertex_descriptor cluster_source = boost::add_vertex(graph);
        Vertex_descriptor cluster_sink = boost::add_vertex(graph);
        std::vector<Vertex_descriptor> inserted_vertices;
        // For E-Data
        // add every facet as a vertex to graph, put edges to source-sink vertices
        for(std::size_t vertex_i = 0; vertex_i <  probability_matrix[0].size();
            ++vertex_i) {
          Vertex_descriptor new_vertex = boost::add_vertex(graph);
          inserted_vertices.push_back(new_vertex);

          double source_weight = probability_matrix[alpha][vertex_i];
          double sink_weight = (labels[vertex_i] == alpha) ?
                               (std::numeric_limits<double>::max)() :
                               probability_matrix[labels[vertex_i]][vertex_i];

          add_edge_and_reverse(cluster_source, new_vertex, source_weight, 0.0, graph);
          add_edge_and_reverse(new_vertex, cluster_sink, sink_weight, 0.0, graph);
        }
        // For E-Smooth
        // add edge between every vertex,
        std::vector<double>::const_iterator weight_it = edge_weights.begin();
        for(std::vector<std::pair<int, int> >::const_iterator edge_it = edges.begin();
            edge_it != edges.end();
            ++edge_it, ++weight_it) {
          Vertex_descriptor v1 = inserted_vertices[edge_it->first];
          Vertex_descriptor v2 = inserted_vertices[edge_it->second];
          add_edge_and_reverse(v1, v2, *weight_it, *weight_it, graph);
        }

        double flow = boost::boykov_kolmogorov_max_flow(graph, cluster_source,
                      cluster_sink);
        if(min_cut - flow < flow * 1e-3) {
          continue;
        }
        std::cout << "prev flow: " << min_cut << "new flow: " << flow << std::endl;
        min_cut = flow;
        success = true;
        //update labeling
        for(std::size_t vertex_i = 0; vertex_i < inserted_vertices.size(); ++vertex_i) {
          boost::default_color_type color = boost::get(boost::vertex_color, graph,
                                            inserted_vertices[vertex_i]);
          if(labels[vertex_i] != alpha && color == ColorTraits::white()) { //new comers
            labels[vertex_i] = alpha;
          }
        }

      }
    } while(success);
  }

  void apply_alpha_expansion_2(const std::vector<std::pair<int, int> >& edges,
                               const std::vector<double>& edge_weights,
                               const std::vector<std::vector<double> >& probability_matrix,
                               std::vector<int>& labels) {
    int number_of_clusters = probability_matrix.size();
    double min_cut = (std::numeric_limits<double>::max)();
    bool success;
    do {
      success = false;
      for(int alpha = 0; alpha < number_of_clusters; ++alpha) {
        Graph graph;
        Vertex_descriptor cluster_source = boost::add_vertex(graph);
        Vertex_descriptor cluster_sink = boost::add_vertex(graph);
        std::vector<Vertex_descriptor> inserted_vertices;
        // For E-Data
        // add every facet as a vertex to graph, put edges to source-sink vertices
        for(std::size_t vertex_i = 0; vertex_i <  probability_matrix[0].size();
            ++vertex_i) {
          Vertex_descriptor new_vertex = boost::add_vertex(graph);
          inserted_vertices.push_back(new_vertex);

          double source_weight = probability_matrix[alpha][vertex_i];
          double sink_weight = (labels[vertex_i] == alpha) ?
                               (std::numeric_limits<double>::max)() :
                               probability_matrix[labels[vertex_i]][vertex_i];

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
            double w1 = label_1 == alpha ? 0 : *weight_it;
            add_edge_and_reverse(v1, v2, w1, w1, graph);
          } else {
            Vertex_descriptor inbetween = boost::add_vertex(graph);
            add_edge_and_reverse(inbetween, cluster_sink, *weight_it, 0.0, graph);

            double w1 = label_1 == alpha ? 0 : *weight_it;
            double w2 = label_2 == alpha ? 0 : *weight_it;

            add_edge_and_reverse(inbetween, v1, w1, w1, graph);
            add_edge_and_reverse(inbetween, v2, w2, w2, graph);
          }
        }

        double flow = boost::boykov_kolmogorov_max_flow(graph, cluster_source,
                      cluster_sink);
        if(min_cut - flow < flow * 1e-3) {
          continue;
        }
        std::cout << "prev flow: " << min_cut << "new flow: " << flow << std::endl;
        min_cut = flow;
        success = true;
        //update labeling
        for(std::size_t vertex_i = 0; vertex_i < inserted_vertices.size(); ++vertex_i) {
          boost::default_color_type color = boost::get(boost::vertex_color, graph,
                                            inserted_vertices[vertex_i]);
          if(labels[vertex_i] != alpha && color == ColorTraits::white()) { //new comers
            labels[vertex_i] = alpha;
          }
        }

      }
    } while(success);
  }
  void apply_simple_cut(const std::vector<std::pair<int, int> >& edges,
                        const std::vector<double>& edge_weights, std::vector<int>& labels,
                        const std::vector<std::vector<double> >& probability_matrix,
                        std::vector<int>& center_ids) {
    //Graph graph(edges.begin(), edges.end(), vertex_weights.size(), edge_weights.size());
    Graph graph;
    Vertex_descriptor cluster_source = boost::add_vertex(graph);
    Vertex_descriptor cluster_sink = boost::add_vertex(graph);
    std::vector<Vertex_descriptor> inserted_vertices;
    //inserted_vertices.reserve(vertex_weights.size());
    for(std::size_t vertex_i = 0; vertex_i <  probability_matrix[0].size();
        ++vertex_i) {
      Vertex_descriptor new_vertex = boost::add_vertex(graph);
      inserted_vertices.push_back(new_vertex);

      Edge_descriptor source_e, source_e_rev, sink_e, sink_e_rev;
      bool source_e_added, source_e_rev_added, sink_e_added, sink_e_rev_added;

      tie(source_e, source_e_added) = boost::add_edge(cluster_source, new_vertex,
                                      graph);
      tie(source_e_rev, source_e_rev_added) = boost::add_edge(new_vertex,
                                              cluster_source, graph);

      tie(sink_e, sink_e_added) = boost::add_edge(new_vertex, cluster_sink, graph);
      tie(sink_e_rev, sink_e_rev_added) = boost::add_edge(cluster_sink, new_vertex,
                                          graph);

      CGAL_assertion(source_e_added && source_e_rev_added &&  sink_e_added
                     && sink_e_rev_added);

      //put vertex capacities (to edges between itself and source & sink)
      boost::put(boost::edge_capacity, graph, source_e,
                 probability_matrix[0][vertex_i]);
      boost::put(boost::edge_capacity, graph, source_e_rev, 0);
      boost::put(boost::edge_capacity, graph, sink_e,
                 probability_matrix[1][vertex_i]);
      boost::put(boost::edge_capacity, graph, sink_e_rev, 0);
      //map reverse edges
      boost::put(boost::edge_reverse, graph, source_e, source_e_rev);
      boost::put(boost::edge_reverse, graph, source_e_rev, source_e);
      boost::put(boost::edge_reverse, graph, sink_e, sink_e_rev);
      boost::put(boost::edge_reverse, graph, sink_e_rev, sink_e);
    }
    std::vector<double>::const_iterator weight_it = edge_weights.begin();
    for(std::vector<std::pair<int, int> >::const_iterator edge_it = edges.begin();
        edge_it != edges.end();
        ++edge_it, ++weight_it) {
      Vertex_descriptor v1 = inserted_vertices[edge_it->first];
      Vertex_descriptor v2 = inserted_vertices[edge_it->second];

      bool edge_added, edge_rev_added;
      Edge_descriptor edge, edge_rev;
      tie(edge, edge_added) = boost::add_edge(v1, v2, graph);
      tie(edge_rev, edge_rev_added) = boost::add_edge(v2, v1, graph);

      CGAL_assertion(edge_added && edge_rev_added);

      //put edge capacities
      boost::put(boost::edge_capacity, graph, edge,  *weight_it);
      boost::put(boost::edge_capacity, graph, edge_rev,  *weight_it);

      //map reverse edges
      boost::put(boost::edge_reverse, graph, edge, edge_rev);
      boost::put(boost::edge_reverse, graph, edge_rev, edge);
    }

    double flow = boost::boykov_kolmogorov_max_flow(graph, cluster_source,
                  cluster_sink);

    for(std::vector<Vertex_descriptor>::iterator vertex_it =
          inserted_vertices.begin();
        vertex_it != inserted_vertices.end(); ++vertex_it) {
      boost::default_color_type color = boost::get(boost::vertex_color, graph,
                                        *vertex_it);
      if(color == ColorTraits::black()) { // in sink
        center_ids.push_back(1);
      } else {
        center_ids.push_back(0);
      }
    }
  }
};
}//namespace internal
}//namespace CGAL
#endif //CGAL_ALPHA_EXPANSION_GRAPH_CUT_H