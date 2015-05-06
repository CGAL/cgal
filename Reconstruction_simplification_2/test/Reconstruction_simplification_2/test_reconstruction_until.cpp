// test_reconstruction_until.cpp

//----------------------------------------------------------
// Test the cgal environment for Reconstruction_simplification_2
//----------------------------------------------------------

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Reconstruction_simplification_2.h>

#include <fstream>

#include <iostream>
#include <iterator>
#include <list>
#include <cassert>

#include "testing_tools.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2                                          Point;
typedef K::FT                                         		FT;
typedef K::Segment_2 										Segment;


int main ()
{

    std::list<Point> points;

    //use the stair example for testing
	load_xy_file_points<Point>("data/stair-noise00.xy", points);

    CGAL::Reconstruction_simplification_2<K> rs2(points);

    rs2.run_until(9);

    rs2.print_stats_debug();

    std::vector<Point> isolated_points;
	std::vector<Segment> edges;

	rs2.extract_list_output(std::back_inserter(isolated_points), std::back_inserter(edges));

    std::cout << "isolated_points " << isolated_points.size() << std::endl;
    std::cout << "edges " << edges.size() << std::endl;

    assert(isolated_points.size() == 0);
    assert(edges.size() == 8);
}
