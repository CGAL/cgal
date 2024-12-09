/** @file yao_exact.cpp
 *
 * A test application that constructs a Yao graph exactly with the vertex list given
 * in a file, and then generates the Gnuplot files to plot the constructed Yao graph.
 *
 * Authors: Weisheng Si, Quincy Tse, Western Sydney University
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_root_of.h>
#include <CGAL/Construct_yao_graph_2.h>
#include <CGAL/gnuplot_output_2.h>

// select the kernel type
typedef CGAL::Exact_predicates_exact_constructions_kernel_with_root_of   Kernel;
typedef Kernel::Point_2                   Point_2;
typedef Kernel::Direction_2               Direction_2;
typedef boost::adjacency_list<boost::setS,
                              boost::vecS,
                              boost::undirectedS,
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
  CGAL::Construct_yao_graph_2<Kernel, Graph> yao(k, initial_direction);
  // create an adjacency_list object
  Graph g;
  // construct the yao graph on the vertex list
  yao(input_begin, input_end, g);

  // obtain the number of vertices in the constructed graph
  boost::graph_traits<Graph>::vertices_size_type n = boost::num_vertices(g);

  // generate gnuplot files for plotting this graph
  std::string fileprefix = "y" + std::to_string(k) + "n" + std::to_string(n);
  CGAL::gnuplot_output_2(g, fileprefix);

  return 0;
}
