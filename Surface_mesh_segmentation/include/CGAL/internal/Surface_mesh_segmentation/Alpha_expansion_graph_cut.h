#ifndef CGAL_SURFACE_MESH_SEGMENTATION_ALPHA_EXPANSION_GRAPH_CUT_H
#define CGAL_SURFACE_MESH_SEGMENTATION_ALPHA_EXPANSION_GRAPH_CUT_H

#include <CGAL/assertions.h>
#include <CGAL/Timer.h>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#include <boost/math/distributions/normal.hpp>

#include <iostream>
#include <fstream>
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

  Alpha_expansion_graph_cut(
    std::vector<std::pair<int, int> >& edges,
    const std::vector<double>& edge_weights,
    const std::vector<std::vector<double> >& probability_matrix,
    std::vector<int>& labels, double* result = NULL) {
    double min_cut = apply_alpha_expansion(edges, edge_weights, probability_matrix,
                                           labels);
    if(result != NULL) {
      *result = min_cut;
    }
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

  double apply_alpha_expansion(const std::vector<std::pair<int, int> >& edges,
                               const std::vector<double>& edge_weights,
                               const std::vector<std::vector<double> >& probability_matrix,
                               std::vector<int>& labels) {
    std::ofstream log_file("log_file.txt");

    // logging input
    log_file << "edges: " << std::endl;
    for(std::vector<std::pair<int, int> >::const_iterator edge_it = edges.begin();
        edge_it != edges.end();
        ++edge_it) {
      log_file << edge_it->first << " " << edge_it->second << std::endl;
    }
    log_file << "edge weights: " << std::endl;
    for(std::vector<double>::const_iterator w_it = edge_weights.begin();
        w_it != edge_weights.end(); ++w_it) {
      log_file << (*w_it) << std::endl;
    }
    log_file << "prob matrix: " << std::endl;
    for(std::vector< std::vector<double> >::const_iterator v_it =
          probability_matrix.begin();
        v_it != probability_matrix.end(); ++v_it) {
      for(std::vector<double>::const_iterator p_it = v_it->begin();
          p_it != v_it->end(); ++p_it) {
        log_file << (*p_it) << " ";
      }
      log_file << std::endl;
    }
    log_file << "labels-input:" << std::endl;
    for(std::vector<int>::const_iterator l_it = labels.begin();
        l_it != labels.end(); ++l_it) {
      log_file << (*l_it) << std::endl;
    }

    ////////////////////////////////////////////////////////////
    int number_of_clusters = probability_matrix.size();
    double min_cut = (std::numeric_limits<double>::max)();
    bool success;
    Timer gt;
    gt.start();
    do {
      success = false;
      for(int alpha = 0; alpha < number_of_clusters; ++alpha) {
        Timer t;
        t.start();
        Graph graph;
#if 0
        graph.m_vertices.reserve(edges.size() + labels.size()); //not documented!
        // in order to see effect of pre-allocation of vector with maximum size
#endif
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
          // since it is expansion move, current alpha labeled vertices will be assigned to alpha again,
          // making sink_weight 'infinity' guarantee this.
          double sink_weight = (labels[vertex_i] == alpha) ?
                               (std::numeric_limits<double>::max)() :
                               probability_matrix[labels[vertex_i]][vertex_i];

          add_edge_and_reverse(cluster_source, new_vertex, source_weight, 0.0, graph);
          add_edge_and_reverse(new_vertex, cluster_sink, sink_weight, 0.0, graph);
        }
        //std::cout << "vertex time: " << t.time() << std::endl;
        t.reset();
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
            // actually no need for this, since two alpha labeled vertices will not be seperated
            // (their edges between sink is infitinity)
            double w1 = (label_1 == alpha) ? 0 : *weight_it;
            add_edge_and_reverse(v1, v2, w1, w1, graph);
          } else {
            Vertex_descriptor inbetween = boost::add_vertex(graph);

            double w1 = (label_1 == alpha) ? 0 : *weight_it;
            double w2 = (label_2 == alpha) ? 0 : *weight_it;
            add_edge_and_reverse(inbetween, v1, w1, w1, graph);
            add_edge_and_reverse(inbetween, v2, w2, w2, graph);
            add_edge_and_reverse(inbetween, cluster_sink, *weight_it, 0.0, graph);
          }
        }
        //std::cout << "edge time: " << t.time() << std::endl;
        t.reset();
        double flow = boost::boykov_kolmogorov_max_flow(graph, cluster_source,
                      cluster_sink);
        //std::cout << "flow time: " << t.time() << std::endl;
        if(min_cut - flow < flow * 1e-10) {
          continue;
        }

        log_file << "prev flow: " << min_cut << " new flow: " << flow << std::endl;
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
    log_file << "labels-output:" << std::endl;
    for(std::vector<int>::const_iterator l_it = labels.begin();
        l_it != labels.end(); ++l_it) {
      log_file << (*l_it) << std::endl;
    }

    std::cout << "Graph-cut time: " << gt.time() <<  std::endl;
    return min_cut;
  }

  double apply_alpha_expansion_2(const std::vector<std::pair<int, int> >& edges,
                                 const std::vector<double>& edge_weights,
                                 const std::vector<std::vector<double> >& probability_matrix,
                                 std::vector<int>& labels) {
    int number_of_clusters = probability_matrix.size();
    double min_cut = (std::numeric_limits<double>::max)();
    bool success;
    Timer gt;
    gt.start();
    double total_time = 0.0;
    do {
      success = false;
      for(int alpha = 0; alpha < number_of_clusters; ++alpha) {
        Timer t;
        t.start();
        Graph graph(edges.size() + labels.size() +
                    2); // allocate using maximum possible size.

        Vertex_descriptor cluster_source = 0;
        Vertex_descriptor cluster_sink = 1;

        // For E-Data
        // add every facet as a vertex to graph, put edges to source-sink vertices
        for(std::size_t vertex_i = 0; vertex_i <  probability_matrix[0].size();
            ++vertex_i) {
          Vertex_descriptor new_vertex = vertex_i + 2;

          double source_weight = probability_matrix[alpha][vertex_i];
          // since it is expansion move, current alpha labeled vertices will be assigned to alpha again,
          // making sink_weight 'infinity' guarantee this.
          double sink_weight = (labels[vertex_i] == alpha) ?
                               (std::numeric_limits<double>::max)() :
                               probability_matrix[labels[vertex_i]][vertex_i];

          add_edge_and_reverse(cluster_source, new_vertex, source_weight, 0.0, graph);
          add_edge_and_reverse(new_vertex, cluster_sink, sink_weight, 0.0, graph);
        }
        total_time += t.time();
        //std::cout << "vertex time: " << t.time() << std::endl;
        t.reset();
        // For E-Smooth
        // add edge between every vertex,
        int b_index = labels.size() + 2;
        std::vector<double>::const_iterator weight_it = edge_weights.begin();
        for(std::vector<std::pair<int, int> >::const_iterator edge_it = edges.begin();
            edge_it != edges.end();
            ++edge_it, ++weight_it) {
          Vertex_descriptor v1 = edge_it->first + 2, v2 = edge_it->second + 2;
          int label_1 = labels[edge_it->first], label_2 = labels[edge_it->second];
          if(label_1 == label_2) {
            // actually no need for this, since two alpha labeled vertices will not be seperated
            // (their edges between sink is infitinity)
            double w1 = (label_1 == alpha) ? 0 : *weight_it;
            add_edge_and_reverse(v1, v2, w1, w1, graph);
          } else {
            Vertex_descriptor inbetween = b_index++;

            double w1 = (label_1 == alpha) ? 0 : *weight_it;
            double w2 = (label_2 == alpha) ? 0 : *weight_it;
            add_edge_and_reverse(inbetween, v1, w1, w1, graph);
            add_edge_and_reverse(inbetween, v2, w2, w2, graph);
            add_edge_and_reverse(inbetween, cluster_sink, *weight_it, 0.0, graph);
          }
        }
        total_time += t.time();
        //std::cout << "edge time: " << t.time() << std::endl;
        t.reset();

        double flow = boost::boykov_kolmogorov_max_flow(graph, cluster_source,
                      cluster_sink);

        total_time += t.time();
        //std::cout << "flow time: " << t.time() << std::endl;
        t.reset();
        if(min_cut - flow < flow * 1e-10) {
          continue;
        }
        //std::cout << "prev flow: " << min_cut << " new flow: " << flow << std::endl;
        min_cut = flow;
        success = true;
        //update labeling

        for(std::size_t vertex_i = 0; vertex_i < labels.size(); ++vertex_i) {
          Vertex_descriptor v = vertex_i + 2;
          boost::default_color_type color = boost::get(boost::vertex_color, graph, v);
          if(labels[vertex_i] != alpha
              && color == ColorTraits::white()) { //new comers (expansion occurs)
            labels[vertex_i] = alpha;
          }
        }
        //std::cout << "label time: " << t.time() << std::endl;
        t.reset();
      }
    } while(success);
    std::cout << "Graph-cut time: " << gt.time() << std::endl;
    return min_cut;
  }

  double apply_alpha_expansion_3(std::vector<std::pair<int, int> >& edges,
                                 const std::vector<double>& edge_weights,
                                 const std::vector<std::vector<double> >& probability_matrix,
                                 std::vector<int>& labels) {
    int number_of_clusters = probability_matrix.size();
    double min_cut = (std::numeric_limits<double>::max)();
    bool success;
    Timer gt;
    gt.start();
    int cluster_source = labels.size();
    int cluster_sink = labels.size() + 1;
    for(std::size_t vertex_i = 0; vertex_i <  probability_matrix[0].size();
        ++vertex_i) {
      edges.push_back(std::pair<int, int>(vertex_i, cluster_source));
      edges.push_back(std::pair<int, int>(vertex_i, cluster_sink));
    }
    int static_edge_count = edges.size();
    do {
      success = false;
      for(int alpha = 0; alpha < number_of_clusters; ++alpha) {
        Timer t;
        t.start();
        Graph graph(edges.begin(), edges.end(), labels.size()+2, edges.size());
        std::cout << "vertex time: " << t.time() << std::endl;




      }
    } while(success);
    std::cout << "tot time: " << gt.time() <<  std::endl;
    return min_cut;
  }
};
}//namespace internal
}//namespace CGAL
#endif //CGAL_SURFACE_MESH_SEGMENTATION_ALPHA_EXPANSION_GRAPH_CUT_H