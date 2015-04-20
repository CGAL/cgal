// test_basic.cpp

//----------------------------------------------------------
// Test the cgal environment for Reconstruction_simplification_2
//----------------------------------------------------------

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Reconstruction_simplification_2.h>
#include "testing_tools.h"

#include<list>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2                                          Point;
typedef K::FT                                         		FT;



int main ()
{
    std::list<Point> points;
	//use the stair example for testing
	load_xy_file_points<Point>("data/stair-noise00.xy", points);

    CGAL::Reconstruction_simplification_2<K> rs2(points.begin(), points.end());

    rs2.run(100); //100 steps

    rs2.print_stats_debug();
}
