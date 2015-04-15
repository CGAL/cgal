// reconstruction_simplification_2_example.cpp


//----------------------------------------------------------
// Simple example for Reconstruction_simplification_2
//----------------------------------------------------------

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Reconstruction_simplification_2.h>


#include<fstream>
#include<iostream>
#include <string>
#include <iterator>
#include <utility>      // std::pair

#include <CGAL/property_map.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2                     	                    Point;
typedef K::Segment_2                 						Segment;

typedef K::FT                                         		FT;

typedef std::pair<Point, FT> PointMassPair;
typedef std::list<PointMassPair> PointMassList;

typedef CGAL::First_of_pair_property_map <PointMassPair>  Point_property_map;
typedef CGAL::Second_of_pair_property_map <PointMassPair> Mass_property_map;

typedef CGAL::Reconstruction_simplification_2<K, Point_property_map, Mass_property_map> Rs_2;




void load_xy_file(const std::string& fileName, PointMassList& points)
{
   std::ifstream ifs(fileName);
   Point point;
   unsigned int nb = 0;
   while (ifs >> point)
   {
	   points.push_back(std::make_pair(point, 1));
   }
   ifs.close();
}



int main ()
{

	PointMassList points;

	load_xy_file("data/stair-noise00.xy", points);

	Point_property_map point_pmap;
	Mass_property_map  mass_pmap;

	Rs_2 rs2(points.begin(), points.end(), point_pmap, mass_pmap);

	rs2.run(100); //100 steps

 	std::vector<Point> isolated_points;
	std::vector<Segment> edges;


	rs2.extract_list_output(std::back_inserter(isolated_points), std::back_inserter(edges));

    std::cerr << "Isolated Vertices" << std::endl;
    for (std::vector<Point>::iterator it = isolated_points.begin();
			it != isolated_points.end(); it++) {
		std::cout  <<  *it << std::endl;
	}

	std::cerr << "Edges" << std::endl;
	for (std::vector<Segment>::iterator it = edges.begin();
			it != edges.end(); it++) {
  		std::cout << *it << std::endl;
     }

}
