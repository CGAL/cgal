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

typedef CGAL::Reconstruction_simplification_2<K> Rs_2;




void load_xy_file(const std::string& fileName, std::list<Point>& points)
{
   std::ifstream ifs(fileName);
   Point point;
   unsigned int nb = 0;
   while (ifs >> point)
   {
	   points.push_back(point);
   }
   ifs.close();
}



int main ()
{

    std::list<Point> points;

	load_xy_file("data/stair-noise00.xy", points);

	Rs_2 rs2(points.begin(), points.end());

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
