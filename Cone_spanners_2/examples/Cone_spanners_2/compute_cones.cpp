/** @file compute_cones.cpp
 * An example application that computes cones given the number of cones (k)
 * and the initial direction. If exact computation is used, for any k<=28,
 * the computation can be done successfully; for any k>28, the computation cannot be completed
 * because CORE::Expr exceeds its limit. We don't experiment with LEDA::real.
 * We believe k<=28 suffices most applications. Also, if inexact computation
 * is used, the computation will be successful for any k>1.
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
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Compute_cone_boundaries_2.h>

// select the kernel type
typedef CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt   Kernel;
typedef Kernel::Point_2                   Point_2;
typedef Kernel::Direction_2               Direction_2;

int main(int argc, char ** argv) {

    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " <no. of cones> [<direction-x> <direction-y>]" << std::endl;
        return 1;
    }

    unsigned long k = atol(argv[1]);
    if (k<2) {
        std::cout << "The number of cones should be larger than 1!" << std::endl;
        return 1;
    }

    Direction_2 initial_direction;
    if (argc == 2)
        initial_direction = Direction_2(1, 0);  // default initial_direction
    else if (argc == 4)
        initial_direction = Direction_2(atof(argv[2]), atof(argv[3]));
    else {
        std::cout << "Usage: " << argv[0] << " <no. of cones> [<direction-x> <direction-y>]" << std::endl;
        return 1;
    }

    // construct the functor, with no parameters needed for the constructor
	CGAL::Compute_cone_boundaries_2<Kernel> cones;
	// create the vector rays to store the results; it should contain no elements initially.
	std::vector<Direction_2> rays;
	// compute the cone boundaries and store them in rays
	cones(k, initial_direction, rays);

	// display the computed rays
    for (int i=0; i<k; i++) 
		std::cout << "Ray " << i << ": " << rays[i] << std::endl;

    return 0;
}
