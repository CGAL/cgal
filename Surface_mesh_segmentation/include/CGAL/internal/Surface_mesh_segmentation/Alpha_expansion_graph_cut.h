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

  Alpha_expansion_graph_cut(std::vector<std::pair<int, int> >& edges,
                            std::vector<double>& edge_weights, std::vector<int>& labels,
                            std::vector<std::vector<double> >& probability_matrix,
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
    std::vector<double>::iterator weight_it = edge_weights.begin();
    for(std::vector<std::pair<int, int> >::iterator edge_it = edges.begin();
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
      if(labels[edge_it->first] == labels[edge_it->second]) {
        boost::put(boost::edge_capacity, graph, edge,  *weight_it);
        boost::put(boost::edge_capacity, graph, edge_rev,  *weight_it);
      } else {
        boost::put(boost::edge_capacity, graph, edge, *weight_it);
        boost::put(boost::edge_capacity, graph, edge_rev, *weight_it);
      }

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
      if(color == ColorTraits::black() ) { // in sink
        center_ids.push_back(1);
      } else {
        center_ids.push_back(0);
      }
    }
    std::cout << flow << std::endl;
    //for(std::pair<Vertex_iterator, Vertex_iterator> pair_v = boost::vertices(graph);
    //    pair_v.first != pair_v.second; ++(pair_v.first))
    //{
    //    if(cluster_source == (*pair_v.first) || cluster_sink == (*pair_v.first)) { continue; }
    //
    //}
  }
};
}//namespace internal
}//namespace CGAL
#endif //CGAL_ALPHA_EXPANSION_GRAPH_CUT_H