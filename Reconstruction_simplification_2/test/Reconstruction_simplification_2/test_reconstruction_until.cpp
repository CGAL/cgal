// test_reconstruction_until.cpp

//----------------------------------------------------------
// Test the cgal environment for Reconstruction_simplification_2
//----------------------------------------------------------

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Reconstruction_simplification_2.h>
#include <CGAL/List_output.h>

#include <fstream>

#include<iostream>
#include <string>
#include <iterator>
#include <utility>      // std::pair
#include <cassert>

#include <CGAL/property_map.h>
#include <CGAL/value_type_traits.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2                                          Point;
typedef K::FT                                         		FT;
typedef K::Segment_2 										Segment;


typedef std::pair<Point, FT> PointMassPair;
typedef std::list<PointMassPair> PointMassList;
typedef PointMassList::const_iterator InputIterator;
typedef CGAL::value_type_traits<InputIterator>::type MassPoint;
typedef CGAL::First_of_pair_property_map <PointMassPair> PointPMap;
typedef CGAL::Second_of_pair_property_map <PointMassPair> MassPMap;


PointMassList* load_xy_file(const std::string& fileName);
PointMassList* simple_point_set();

int main ()
{

	//use the stair example for testing
	PointMassList points = *(load_xy_file("data/stair-noise00.xy"));

    PointPMap point_pmap;
    MassPMap  mass_pmap;

    MassPoint mp;

    CGAL::Reconstruction_simplification_2<K, InputIterator, PointPMap, MassPMap>
    	rs2(points.begin(), points.end(), point_pmap, mass_pmap);

    rs2.initialize();

    rs2.reconstruct_until(9);

    rs2.print_stats_debug();

    std::vector<Point> isolated_points;
	std::vector<Segment> edges;

	typedef std::back_insert_iterator<std::vector<Point> >   Point_it;
	typedef std::back_insert_iterator<std::vector<Segment> > Edge_it;

	Point_it point_it(isolated_points);
	Edge_it  edge_it(edges);

	CGAL::List_output<K, Point_it, Edge_it> list_output(point_it, edge_it);
    rs2.extract_solid_elements(list_output);


    std::cout << "isolated_points " << isolated_points.size() << std::endl;
    std::cout << "edges " << edges.size() << std::endl;

    assert(isolated_points.size() == 0);
    assert(edges.size() == 8);
}


PointMassList* simple_point_set() {

	 PointMassList *points = new PointMassList();

	 points->push_back(std::make_pair(Point(0.0001, 0.00012), 1));
	 points->push_back(std::make_pair(Point(1.00013, 0.00031), 1));
	 points->push_back(std::make_pair(Point(0.50001, 0.500001), 1));
	 points->push_back(std::make_pair(Point(0.5000011, 1.000012), 1));
	 points->push_back(std::make_pair(Point(0.000012, 2.000001), 1));
	 points->push_back(std::make_pair(Point(1.0000123, 2.0000012), 1));
	 points->push_back(std::make_pair(Point(-0.50000182, 0.50000162), 1));
	 points->push_back(std::make_pair(Point(-0.500001002, 1.000006012), 1));
	 points->push_back(std::make_pair(Point(-1.0000331, 2.00001), 1));
	 points->push_back(std::make_pair(Point(-1.000012, 0.000200012), 1));

    return points;
}


PointMassList* load_xy_file(const std::string& fileName)
{
	PointMassList *points = new PointMassList();
       std::ifstream ifs(fileName);
       std::cerr << "read xy...";
       Point point;
       unsigned int nb = 0;
       while (ifs >> point)
       {
    	   points->push_back(std::make_pair(point, 1));
       }
       std::cerr << "done (" << nb << " points)" << std::endl;
       ifs.close();

       return points;

}
