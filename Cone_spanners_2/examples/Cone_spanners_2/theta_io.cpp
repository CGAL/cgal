/** @file theta_io.cpp
 *
 * An example application that constructs a Theta graph with an input vertex list,
 * and then generates the Gnuplot files to plot the Theta graph.
 */

// authors: Weisheng Si, Quincy Tse
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>
#include <algorithm>

#include <boost/graph/adjacency_list.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Construct_theta_graph_2.h>
#include <CGAL/gnuplot_output_2.h>

// select the kernel type
typedef CGAL::Exact_predicates_inexact_constructions_kernel   Kernel;
typedef Kernel::Point_2                   Point_2;
typedef Kernel::Direction_2               Direction_2;
// it is important that the edgelist is 'setS', such that duplicate edges will be automatically removed.
typedef boost::adjacency_list<boost::setS,
                              boost::vecS,
                              boost::undirectedS,
                              Point_2
                             > Graph;

int main(int argc, char ** argv) {

    if (argc < 3) {
        std::cout << "Usage: " << argv[0] << " <no. of cones> <input filename> [<direction-x> <direction-y>]" << std::endl;
        return 1;
    }

    unsigned long k = atol(argv[1]);
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
    std::string fileprefix = "t" + std::to_string(k) + "n" + std::to_string(n);
    CGAL::gnuplot_output_2(g, fileprefix);

    return 0;
}
