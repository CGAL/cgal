// reconstruction_simplification_2_output_example.cpp


//----------------------------------------------------------
// Simple output example for Reconstruction_simplification_2
//----------------------------------------------------------


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Reconstruction_simplification_2.h>

#include<fstream>
#include<iostream>
#include <string>
#include <iterator>
#include <utility>      // std::pair


#include <CGAL/property_map.h>
#include <CGAL/value_type_traits.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2                     	                    Point;
typedef K::Segment_2                 						Segment;

typedef K::FT                                         		FT;

typedef std::pair<Point, FT> PointMassPair;
typedef std::list<PointMassPair> PointMassList;

typedef CGAL::First_of_pair_property_map <PointMassPair>  PointPMap;
typedef CGAL::Second_of_pair_property_map <PointMassPair> MassPMap;


typedef CGAL::Reconstruction_simplification_2<K, PointPMap, MassPMap> Rs_2;

typedef Rs_2::Vertex Vertex;

typedef Rs_2::Reconstruction_edge_2 R_edge_2;

typedef CGAL::Reconstruction_triangulation_2<K> Rt_2;

typedef Rt_2::Finite_edges_iterator Finite_edges_iterator;
typedef Rt_2::Vertex_iterator Vertex_iterator;

typedef Rt_2::Edge Edge;


void list_output(Rs_2& rs2);
void tds_output(Rs_2& rs2);
void index_output(Rs_2& rs2);


void print_edge(Edge edge) {
	int i = edge.second;
	const Point& a = edge.first->vertex((i+1)%3)->point();
	const Point& b = edge.first->vertex((i+2)%3)->point();
	std::cout << a << " , " << b  << std::endl;

}

void load_xy_file(const std::string& filename, PointMassList& points)
{
   std::ifstream ifs(filename);
   Point point;
   while (ifs >> point)
	   points.push_back(std::make_pair(point, 1));

   ifs.close();
}

int main ()
{

	PointMassList points;
	load_xy_file("data/stair-noise00.xy", points);

	PointPMap point_pmap;
	MassPMap  mass_pmap;

	Rs_2 rs2(points.begin(), points.end(), point_pmap, mass_pmap);

	rs2.reconstruct(100); //100 steps

	list_output(rs2);
	tds_output(rs2);
	index_output(rs2);
}

void list_output(Rs_2& rs2) {

	std::cout << "(-------------List output---------- )" << std::endl;

	std::vector<Point> isolated_points;
	std::vector<Segment> edges;

	rs2.extract_list_output(std::back_inserter(isolated_points), std::back_inserter(edges));

	for (std::vector<Point>::iterator it = isolated_points.begin();
			it != isolated_points.end(); it++) {
		std::cout  <<  *it << std::endl;
	}

	for (std::vector<Segment>::iterator it = edges.begin();
			it != edges.end(); it++) {
  		std::cout << *it << std::endl;
    }
}

void tds_output(Rs_2& rs2) {

	std::cout << "(-------------Tds output---------- )" << std::endl;

	Rt_2 rt2;
	rs2.extract_tds_output(rt2);

	for (Vertex_iterator vi = rt2.vertices_begin();
					  vi != rt2.vertices_end(); ++vi) {

		FT relevance = (*vi).get_relevance();
		if (relevance > 0) {
			std::cout  <<  *vi << std::endl;
		}
	}

	for (Finite_edges_iterator ei = rt2.finite_edges_begin(); ei != rt2.finite_edges_end(); ++ei) {
		FT relevance = (*ei).first->relevance((*ei).second);
		if (relevance > 0) {
			print_edge(*ei);
		}
	}
}

void index_output(Rs_2& rs2) {

	std::cout << "(-------------Off output---------- )" << std::endl;

	rs2.extract_index_output(std::cout);
}
