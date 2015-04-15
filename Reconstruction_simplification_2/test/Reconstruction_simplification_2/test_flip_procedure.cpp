
// test_flip_procedure.cpp

//----------------------------------------------------------
// Test the cgal environment for Reconstruction_simplification_2
//----------------------------------------------------------

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Reconstruction_simplification_2.h>
#include "testing_tools.h"

#include<iostream>
#include <string>
#include <iterator>
#include <utility>      // std::pair
#include <cassert>

#include <CGAL/property_map.h>



typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2                                          Point;
typedef K::FT                                         		FT;

typedef std::pair<Point, FT> PointMassPair;
typedef std::list<PointMassPair> PointMassList;

typedef CGAL::First_of_pair_property_map <PointMassPair> Point_property_map;
typedef CGAL::Second_of_pair_property_map <PointMassPair> Mass_property_map;


int main ()
{
	PointMassList points;
	//use the stair example for testing
	load_xy_file<PointMassList, Point>("data/stair-noise00.xy", points);

    Point_property_map point_pmap;
    Mass_property_map  mass_pmap;


    for (int i = 1; i <= points.size();) {

		CGAL::Reconstruction_simplification_2<K, Point_property_map, Mass_property_map>
			rs2(points.begin(), points.end(), point_pmap, mass_pmap);

		rs2.run_until(i);

		rs2.print_stats_debug();

		assert(rs2.get_vertex_count() == i);

		i = i + 20;

    }
}
