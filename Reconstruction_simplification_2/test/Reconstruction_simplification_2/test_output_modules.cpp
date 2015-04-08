// test_output_modules.cpp

//----------------------------------------------------------
// Test the cgal environment for Reconstruction_simplification_2
//----------------------------------------------------------

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Reconstruction_simplification_2.h>

#include<fstream>
#include<iostream>
#include <string>
#include <cassert>
#include <iterator>
#include <utility>      // std::pair


#include <CGAL/property_map.h>
#include "testing_tools.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2                     	                    Point;
typedef K::Segment_2                 						Segment;

typedef K::FT                                         		FT;

typedef std::pair<Point, FT> PointMassPair;
typedef std::list<PointMassPair> PointMassList;

typedef CGAL::First_of_pair_property_map <PointMassPair>  Point_property_map;
typedef CGAL::Second_of_pair_property_map <PointMassPair> Mass_property_map;


typedef CGAL::Reconstruction_simplification_2<K, Point_property_map, Mass_property_map> Rs_2;

typedef Rs_2::Vertex Vertex;

typedef Rs_2::Reconstruction_edge_2 R_edge_2;

typedef CGAL::Reconstruction_triangulation_2<K> Rt_2;

typedef Rt_2::Finite_edges_iterator Finite_edges_iterator;
typedef Rt_2::Vertex_iterator Vertex_iterator;

typedef Rt_2::Edge Edge;



void test_list_output(Rs_2& rs2);
void test_index_output(Rs_2& rs2);
void test_tds_output(Rs_2& rs2);



void print_edge(Edge edge) {
	int i = edge.second;
	Point a = edge.first->vertex((i+1)%3)->point();
	Point b = edge.first->vertex((i+2)%3)->point();
	std::cout << a << " , " << b <<  std::endl;
}

int main ()
{

	PointMassList points;
	//use the stair example for testing
	load_xy_file<PointMassList, Point>("data/stair-noise00.xy", points);

	Point_property_map point_pmap;
	Mass_property_map  mass_pmap;

	Rs_2 rs2(points.begin(), points.end(), point_pmap, mass_pmap);

	rs2.reconstruct(100); //100 steps

	test_list_output(rs2);
	test_index_output(rs2);
	test_tds_output(rs2);

}


void test_index_output(Rs_2& rs2) {

	std::cout <<"(-------------Off OUTPUT---------- )" << std::endl;

    rs2.extract_index_output(std::cout);

    //print

    //test cardinalities
    std::ostringstream buffer;
    rs2.extract_index_output(buffer);

    std::stringstream stream(buffer.str());
	std::vector<std::string> res;
	while (1){
		std::string line;
		std::getline(stream,line);
		if (!stream.good())
			break;
		res.push_back(line);
	}

	assert(res.size() == 110);

	assert(res.front() == "OFF 60 0 31");

	for (int i = 61; i < 79; i++) {
		assert(res[i].substr(0,2) == "1 ");
	}

	for (int i = 79; i < 110; i++) {
		assert(res[i].substr(0,2) == "2 ");
	}
}

void test_list_output(Rs_2& rs2) {

	std::cout <<"(-------------List OUTPUT---------- )" << std::endl;

	std::vector<Point> isolated_points;
	std::vector<Segment> edges;

	rs2.extract_list_output(std::back_inserter(isolated_points), std::back_inserter(edges));

    int vertex_count = 0;
	for (std::vector<Point>::iterator it = isolated_points.begin();
			it != isolated_points.end(); it++) {
		vertex_count++;
  		std::cout  <<  *it << std::endl;
	}
	assert(vertex_count == 18);

	int edge_count = 0;
	for (std::vector<Segment>::iterator it = edges.begin();
			it != edges.end(); it++) {
  		std::cout << *it << std::endl;
  		edge_count++;
    }
	assert(edge_count == 31);


}
