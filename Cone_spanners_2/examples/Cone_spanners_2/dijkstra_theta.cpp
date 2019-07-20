#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Construct_theta_graph_2.h>
#include <CGAL/property_map.h>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

// select the kernel type
typedef CGAL::Exact_predicates_inexact_constructions_kernel   Kernel;
typedef Kernel::Point_2                   Point_2;
typedef Kernel::Direction_2               Direction_2;

/* define the struct for edge property */
struct Edge_property {
  /* record the Euclidean length of the edge */
  double euclidean_length;
};

// define the Graph (e.g., to be undirected,
// and to use Edge_property as the edge property, etc.)
typedef boost::adjacency_list<boost::listS,
                              boost::vecS,
                              boost::undirectedS,
                              Point_2,
                              Edge_property>  Graph;

int main(int argc, char ** argv)
{

  if (argc != 3) {
    std::cout << "Usage: " << argv[0] << " <no. of cones> <input filename>" << std::endl;
    return 1;
  }
  unsigned int k = atoi(argv[1]);
  if (k<2) {
    std::cout << "The number of cones should be larger than 1!" << std::endl;
    return 1;
  }
  // open the file containing the vertex list
  std::ifstream inf(argv[2]);
  if (!inf) {
    std::cout << "Cannot open file " << argv[1] << "!" << std::endl;
    return 1;
  }

  // iterators for reading the vertex list file
  std::istream_iterator< Point_2 > input_begin( inf );
  std::istream_iterator< Point_2 > input_end;

  // initialize the functor
  // If the initial direction is omitted, the x-axis will be used
  CGAL::Construct_theta_graph_2<Kernel, Graph> theta(k);
  // create an adjacency_list object
  Graph g;
  // construct the theta graph on the vertex list
  theta(input_begin, input_end, g);

  // select a source vertex for dijkstra's algorithm
  boost::graph_traits<Graph>::vertex_descriptor v0;
  v0 = vertex(0, g);
  std::cout << "The source vertex is: " << g[v0] << std::endl;

  std::cout << "The index of source vertex is: " << v0 << std::endl;

  // calculating edge length in Euclidean distance and store them in the edge property
  boost::graph_traits<Graph>::edge_iterator ei, ei_end;
  for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
    boost::graph_traits<Graph>::edge_descriptor e = *ei;
    boost::graph_traits<Graph>::vertex_descriptor  u = source(e, g);
    boost::graph_traits<Graph>::vertex_descriptor  v = target(e, g);
    const Point_2& pu = g[u];
    const Point_2& pv = g[v];

    double dist = CGAL::sqrt( CGAL::to_double(CGAL::squared_distance(pu,pv)) );
    g[e].euclidean_length = dist;
    std::cout << "Edge (" << g[u] << ", " << g[v] << "): " << dist << std::endl;
  }

  // calculating the distances from v0 to other vertices
  boost::graph_traits<Graph>::vertices_size_type n = num_vertices(g);
  // vector for storing the results
  std::vector<double> distances(n);
  // Calling the Dijkstra's algorithm implementation from boost.
  boost::dijkstra_shortest_paths(g,
                                 v0,
                                 boost::weight_map(get(&Edge_property::euclidean_length, g)).
                                 distance_map(CGAL::make_property_map(distances)) );

  std::cout << "distances are:" << std::endl;
  for (unsigned int i=0; i < n; ++i) {
    std::cout << "distances[" << i << "] = " << distances[i] << ", (x,y)=" << g[vertex(i, g)];
    std::cout << " at Vertex " << i << std::endl;
  }

  return 0;
}
