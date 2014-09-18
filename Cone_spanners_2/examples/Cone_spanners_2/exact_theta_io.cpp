// Copyright (c) 2013-2014  The University of Western Sydney, Australia.
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
//
// Authors: Weisheng Si, Quincy Tse

/** @file exact_theta_io.cpp
 *
 * An example application that exactly constructs a Theta graph with an input vertex list,
 * and then generates the Gnuplot files to plot the Theta graph.
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>
#include <algorithm>

#include <CGAL/basic.h>

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>

#include "Theta_graph_2.h"
#include "gnuplot_output_2.h"

// select the kernel type
typedef CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt   Kernel;
typedef Kernel::Point_2                   Point_2;
typedef Kernel::Direction_2               Direction_2;
// define the theta graph to use the selected kernel and to be undirected
typedef CGAL::Theta_graph_2<Kernel, boost::undirectedS>      T;
// obtain the graph type by boost::adjacency_list
typedef T::Graph                          Graph;

int main(int argc, char ** argv) {

    if (argc < 3) {
        std::cout << "Usage: " << argv[0] << " <no. of cones> <input filename> [<dx> <dy>]" << std::endl;
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

    Direction_2 startingray;
    if (argc == 3)
        startingray = Direction_2(1, 0);
    else if (argc == 5)
        startingray = Direction_2(atof(argv[3]), atof(argv[4]));
    else {
        std::cout << "Usage: " << argv[0] << " <no. of cones> <input filename> [<dx> <dy>]" << std::endl;
        return 1;
    }

    // iterators for reading the vertex list file
    std::istream_iterator<Point_2> input_begin( inf );
    std::istream_iterator<Point_2> input_end;

    // construct the theta graph on the vertex list
    T t(k, input_begin, input_end, startingray);

    // obtain a reference to the boost::adjacency_list object of the constructed graph
    const Graph& g = t.get_graph();
    // obtain the number of vertices in the constructed graph
    unsigned int n = boost::num_vertices(g);

    // generate gnuplot files for plotting this graph, 'e' stands for 'exact'
    std::string fileprefix = "et" + std::to_string(k) + "n" + std::to_string(n);
    CGAL::gnuplot_output_2(g, fileprefix);

    return 0;
}
