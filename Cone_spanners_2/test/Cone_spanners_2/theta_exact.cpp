/** @file theta_exact.cpp
 * A test application that constructs a Theta graph exactly with the vertex list
 * given in a file, and then generates the Gnuplot files for plotting the constructed Theta graph.
 *
 * Authors: Weisheng Si and Quincy Tse, University of Western Sydney
 */
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>
#include <boost/lexical_cast.hpp>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_root_of.h>
#include <CGAL/Construct_theta_graph_2.h>
#include <CGAL/gnuplot_output_2.h>

// select the kernel type
typedef CGAL::Exact_predicates_exact_constructions_kernel_with_root_of   Kernel;
typedef Kernel::Point_2                   Point_2;
typedef Kernel::Direction_2               Direction_2;
// define the graph type
typedef boost::adjacency_list<boost::listS,
                              boost::vecS,
                              #ifdef CGAL_CXX11
                              boost::undirectedS,
                              #else
                              boost::directedS,
                              #endif
                              Point_2
                             > Graph;

int main(int argc, char ** argv)
{
  if (argc < 3) {
    std::cout << "Usage: " << argv[0] << " <no. of cones> <input filename> [<direction-x> <direction-y>]" << std::endl;
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
    std::cout << "Cannot open file " << argv[2] << "!" << std::endl;
    return 1;
  }

  Direction_2 initial_direction;
  if (argc == 3)
    initial_direction = Direction_2(1, 0);  // default initial_direction
  else if (argc == 5)
    initial_direction = Direction_2(atof(argv[3]), atof(argv[4]));
  else {
    std::cout << "Usage: " << argv[0] << " <no. of cones> <input filename> [<direction-x> <direction-y>]" << std::endl;
    return 1;
  }

  // iterators for reading the vertex list file
  std::istream_iterator<Point_2> input_begin( inf );
  std::istream_iterator<Point_2> input_end;

  // initialize the functor
  CGAL::Construct_theta_graph_2<Kernel, Graph> theta(k, initial_direction);
  // create an adjacency_list object
  Graph g;
  // construct the theta graph on the vertex list
  theta(input_begin, input_end, g);

  // obtain the number of vertices in the constructed graph
  unsigned int n = boost::num_vertices(g);
  // generate gnuplot files for plotting this graph
  std::string file_prefix = "t" + boost::lexical_cast<std::string>(k) + "n" + boost::lexical_cast<std::string>(n);
  CGAL::gnuplot_output_2(g, file_prefix);

  return 0;
}
