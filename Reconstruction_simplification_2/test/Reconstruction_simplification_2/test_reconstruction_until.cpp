// test_reconstruction_until.cpp

//----------------------------------------------------------
// Test the cgal environment for Reconstruction_simplification_2
//----------------------------------------------------------

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Reconstruction_simplification_2.h>

#include <fstream>

#include<iostream>
#include <string>
#include <iterator>
#include <utility>      // std::pair
#include <cassert>

#include <CGAL/property_map.h>
#include "testing_tools.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2                                          Point;
typedef K::FT                                         		FT;
typedef K::Segment_2 										Segment;


typedef std::pair<Point, FT> PointMassPair;
typedef std::list<PointMassPair> PointMassList;

typedef CGAL::First_of_pair_property_map <PointMassPair> PointPMap;
typedef CGAL::Second_of_pair_property_map <PointMassPair> MassPMap;

int main ()
{

	PointMassList points;
	//use the stair example for testing
	load_xy_file<PointMassList, Point>("data/stair-noise00.xy", points);

    PointPMap point_pmap;
    MassPMap  mass_pmap;

    CGAL::Reconstruction_simplification_2<K, PointPMap, MassPMap>
    	rs2(points.begin(), points.end(), point_pmap, mass_pmap);

    rs2.reconstruct_until(9);

    rs2.print_stats_debug();

    std::vector<Point> isolated_points;
	std::vector<Segment> edges;

	rs2.extract_list_output(std::back_inserter(isolated_points), std::back_inserter(edges));

    std::cout << "isolated_points " << isolated_points.size() << std::endl;
    std::cout << "edges " << edges.size() << std::endl;

    assert(isolated_points.size() == 0);
    assert(edges.size() == 8);
}
