// test_list_output.cpp

//----------------------------------------------------------
// Test the cgal environment for Reconstruction_simplification_2
//----------------------------------------------------------

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Reconstruction_simplification_2.h>
#include <CGAL/List_output.h>
#include <CGAL/Off_output.h>
#include <CGAL/Tds_output.h>

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
typedef PointMassList::const_iterator InputIterator;
typedef CGAL::value_type_traits<InputIterator>::type MassPoint;
typedef CGAL::First_of_pair_property_map <PointMassPair>  PointPMap;
typedef CGAL::Second_of_pair_property_map <PointMassPair> MassPMap;
typedef CGAL::Reconstruction_simplification_2 <K, InputIterator,
					PointPMap, MassPMap>::Reconstruction_edge_2 R_edge_2;

typedef CGAL::Reconstruction_simplification_2 <K, InputIterator,
					PointPMap, MassPMap>::Vertex Vertex;

typedef CGAL::Reconstruction_triangulation_2<K> Rt_2;

typedef CGAL::List_output<K>::Output_Vertex_Iterator Output_Vertex_Iterator;
typedef CGAL::List_output<K>::Output_Edge_Iterator   Output_Edge_Iterator;

typedef typename Rt_2::Finite_edges_iterator Finite_edges_iterator;
typedef typename Rt_2::Vertex_iterator Vertex_iterator;



PointMassList* load_xy_file(const std::string& fileName);
PointMassList* simple_point_set();


void print_vertex(Point vertex) {
	std::cout  <<  vertex << std::endl;
}


void print_edge(Segment edge) {
	/*int i = ((edge).edge()).second;
	Point a = ((edge).edge()).first->vertex((i+1)%3)->point();
	Point b = ((edge).edge()).first->vertex((i+2)%3)->point();
	std::cout << "( " << a << " , " << b << " )" << std::endl;*/
	//"( " << (edge).priority()  <<  ")


	std::cout << edge << std::endl;
}

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

    rs2.reconstruct(100); //100 steps

    rs2.print_stats_debug();

    CGAL::List_output<K> list_output;

    rs2.extract_solid_elements(list_output);

    std::cout <<"(-------------List OUTPUT---------- )" << std::endl;


  	for (Output_Vertex_Iterator it = list_output.vertices_start();
			it != list_output.vertices_beyond(); it++) {
  		print_vertex(*it);
   }

	for (Output_Edge_Iterator it = list_output.edges_start();
			it != list_output.edges_beyond(); it++) {
		print_edge(*it);
    }/*

	//-------
	std::cout <<"(-------------OFF OUTPUT----------- )" << std::endl;

    CGAL::Off_output<K> off_output;
    rs2.extract_solid_elements(off_output);
    off_output.get_os_output(std::cout);


	//-------
	std::cout <<"(-------------TRI OUTPUT----------- )" << std::endl;

    CGAL::Tds_output<K> tds_output;
    rs2.extract_solid_elements(tds_output);
    Rt_2 rt2;
    tds_output.extract_reconstruction_tds(rt2);


    for (Vertex_iterator vi = rt2.vertices_begin();
    				  vi != rt2.vertices_end(); ++vi) {

    	FT relevance = (*vi).get_relevance();
		if (relevance <= 0)
			continue;

		print_vertex(*vi);
    }

    for (Finite_edges_iterator ei = rt2.finite_edges_begin(); ei != rt2.finite_edges_end(); ++ei) {
    	FT relevance = (*ei).first->relevance((*ei).second);
    	if (relevance <= 0)
    		continue;

    	std::cout <<  relevance;
    	print_edge(*ei);
    	(*ei).first->relevance((*ei).second);
    }*/
}

PointMassList* simple_point_set() {

	PointMassList *points = new PointMassList();

    points->push_back(std::make_pair(Point(0.1,0.1), 1));
    points->push_back(std::make_pair(Point(0.4,0.1), 1));
    points->push_back(std::make_pair(Point(0.6,0.1), 1));
    points->push_back(std::make_pair(Point(0.9,0.1), 1));
    points->push_back(std::make_pair(Point(0.9,0.4), 1));
    points->push_back(std::make_pair(Point(0.9,0.6), 1));
    points->push_back(std::make_pair(Point(0.9,0.9), 1));
    points->push_back(std::make_pair(Point(0.6,0.9), 1));
    points->push_back(std::make_pair(Point(0.4,0.9), 1));
    points->push_back(std::make_pair(Point(0.1,0.9), 1));
    points->push_back(std::make_pair(Point(0.1,0.6), 1));
    points->push_back(std::make_pair(Point(0.1,0.4), 1));

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
