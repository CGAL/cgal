
// test_flip_procedure.cpp

//----------------------------------------------------------
// Test the cgal environment for Reconstruction_simplification_2
//----------------------------------------------------------

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Reconstruction_simplification_2.h>
#include "testing_tools.h"

#include <iostream>
#include <string>
#include <iterator>
#include <utility>      // std::pair
#include <cassert>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2                                          Point;
typedef K::FT                                         		FT;


int main ()
{
    std::list<Point> points;
	//use the stair example for testing
	load_xy_file_points<Point>("data/stair-noise00.xy", points);


    for (int i = 1; i <= points.size();) {

		CGAL::Reconstruction_simplification_2<K> rs2(points);

		rs2.run_until(i);

		rs2.print_stats_debug();

		assert(rs2.number_of_vertices() == i);

		i = i + 20;

    }
}
