//----------------------------------------------------------
//Pm includes
//----------------------------------------------------------
#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_segment_cached_traits_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/Pm_with_intersections.h>
#include <CGAL/IO/Pm_iostream.h>
#include <CGAL/Pm_trapezoid_ric_point_location.h>
#include <CGAL/Pm_naive_point_location.h>
#include <CGAL/Pm_simple_point_location.h>
#include <CGAL/Pm_walk_along_line_point_location.h>
#include <CGAL/Pm_triangle_point_location.h>
#include <CGAL/Pm_landmarks_point_location.h>
#include <CGAL/Pm_nearest_neighbor.h>
#include "Segment_reader.h"

#include <iostream>
#include <fstream>
#include <stdio.h>

//----------------------------------------------------------
// Pm Types
//----------------------------------------------------------

typedef CGAL::Quotient<CGAL::MP_Float>                  Coord_type;
typedef CGAL::Cartesian<Coord_type>                     Kernel;
typedef CGAL::Arr_segment_cached_traits_2<Kernel>       Traits;

typedef Traits::Point_2                                 Point_2;
typedef Traits::Curve_2                                 Curve_2;

typedef CGAL::Pm_default_dcel<Traits>                   Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>                 Pm;
typedef CGAL::Planar_map_with_intersections_2<Pm>       Pmwx;
typedef CGAL::Pm_trapezoid_ric_point_location<Pm>       Trap_point_location;
typedef CGAL::Pm_naive_point_location<Pm>               Naive_point_location;
typedef CGAL::Pm_simple_point_location<Pm>              Simple_point_location;
typedef CGAL::Pm_walk_along_line_point_location<Pm>     Walk_point_location;
typedef CGAL::Pm_nearest_neighbor<Pm>               Nearest_neighbor;
typedef CGAL::Pm_landmarks_point_location<Pm,Nearest_neighbor>  Landmarks_point_location;
typedef Pm::Face_iterator                               Face_iterator;
typedef Pm::Halfedge_iterator                           Halfedge_iterator;
typedef Pm::Vertex_iterator                             Vertex_iterator;
typedef Pm::Edge_iterator                               Edge_iterator;
typedef Pm::Vertex_handle                               Vertex_handle;
typedef Pm::Halfedge_around_vertex_circulator Halfedge_around_vertex_circulator;
typedef Pm::Ccb_halfedge_circulator                     Ccb_halfedge_circulator;
typedef Pm::Locate_type                                 Locate_type;
typedef std::list<Curve_2>                              Curve_list;
typedef std::list<Point_2>                                Point_list;
typedef Pm::Traits_wrap                  Traits_wrap;
//----------------------------------------------------------

/*!
*/
int main(int argc, char * argv[] )
{
	const char * inp_seg_filename = "temp_input_seg.txt";
	const char * inp_pnt_filename = "temp_input_point.txt";
	const char * out_filename = "temp_results.txt";

	Landmarks_point_location		strategy;
	Pmwx											pm(&strategy);
	Locate_type								lt;
	Halfedge_iterator					e;
	Traits											tr;
	Traits_wrap *							traits = (Traits_wrap*)(&tr);
	Curve_list									curve_list;
	Point_list                                  point_list;

	//read the segments from the first input file
	Segment_reader<Traits> reader;
	CGAL::Bbox_2 bbox(0.0,0.0,0.0,0.0);

	reader.read_data(inp_seg_filename, std::back_inserter(curve_list), bbox);
	std::cout << curve_list.size() << " curves" << std::endl;
	
	//incremental insert
	//Curve_list::const_iterator cv_iter; 
	//for (cv_iter = curve_list.begin(); cv_iter != curve_list.end(); cv_iter++) {
	//	std::cout << "+++ insert cv = " << *cv_iter << std::endl;
	//	pm.insert(*cv_iter) ;
	//}

	//aggregate insert
	pm.insert(curve_list.begin(), curve_list.end());	
	std::cout << "All curves were inserted to the arrangement" << std::endl;
	
	//std::cout << pm;
	//getchar();
	
	 //read points from file into list
 	std::ifstream inp_pnt_file(inp_pnt_filename);
	if (!inp_pnt_file.is_open()) 
	{
		std::cerr << "Cannot open file " << inp_pnt_filename << "!" << std::endl;
		getchar(); return -1;
	}
    
    std::cout << "before read points " << std::endl;
	int points_count;
    inp_pnt_file >> points_count;
    
    for (int i = 0; i < points_count; i++) {
      Coord_type x, y;
      inp_pnt_file >> x >> y;
      Point_2 pnt(x, y);
      point_list.push_back(pnt);
    }

     std::cout << point_list.size() << " points" << std::endl;
	 getchar();
    
    //iterator on all the points got from the input, 
    //go over each one and locate it in the pm .
	Point_list::const_iterator iter;
	int point_index =0;
	for (iter = point_list.begin(); iter != point_list.end(); iter++) 
	{
		Point_2 pnt = *iter;
		std::cout << std::endl << std::endl;
		std::cout << "-------- point number "<< point_index <<" is " << pnt << std::endl << std::endl;
		//if (point_index != 384) {
		//	point_index++;
		//	continue;
		//}
		point_index++;
	
		//locate point
		e = pm.locate(pnt,lt);

		//print output
		if (lt==Pm::UNBOUNDED_FACE) std::cout << "Unbounded face" << std::endl;
		else if (lt==Pm::FACE) std::cout << "Face that is left of :  ";
		else if (lt==Pm::EDGE) std::cout << "EDGE : " ;
		else if (lt==Pm::VERTEX) std::cout << "VERTEX : " ;
		else std::cout << "Unknown locate type" << std::endl;
		//print e
		if (lt != Pm::UNBOUNDED_FACE)
		{
			std::cout << "e = "<< e->source()->point() <<" towards "<< e->target()->point() << std::endl;
			//outfile << e->curve() << std::endl;
			//if (lt == Pm::FACE) {
			//	Ccb_halfedge_circulator curr = e;
			//	curr ++;
			//	while (curr != e) {
			//		outfile << curr->curve() << std::endl;	
			//		curr ++;
			//	}
			//}
		}
	}

	std::cout << "we're done!" << std::endl; 
	getchar();

	return 0;
}

	//house
	//curve_list.push_back(Curve_2(Point_2(0,0), Point_2(4,0)));
	//curve_list.push_back(Curve_2(Point_2(0,0), Point_2(0,4)));
	//curve_list.push_back(Curve_2(Point_2(4,0), Point_2(4,4)));
	//curve_list.push_back(Curve_2(Point_2(0,4), Point_2(4,4)));
	//curve_list.push_back(Curve_2(Point_2(0,4), Point_2(2,8)));
	//curve_list.push_back(Curve_2(Point_2(4,4), Point_2(2,8))); 

	//special case 1
	//curve_list.push_back(Curve_2(Point_2(4,2), Point_2(20,4))); 
	//curve_list.push_back(Curve_2(Point_2(4,2), Point_2(4,4))); 
	//curve_list.push_back(Curve_2(Point_2(20,4), Point_2(4,4))); 
	//curve_list.push_back(Curve_2(Point_2(11,5), Point_2(9,7))); 

	////special case 2
	//curve_list.push_back(Curve_2(Point_2(4,8), Point_2(19,9))); 
	//curve_list.push_back(Curve_2(Point_2(4,10), Point_2(19,9))); 
	//curve_list.push_back(Curve_2(Point_2(4,8), Point_2(4,10))); 
	//curve_list.push_back(Curve_2(Point_2(11,11), Point_2(13,13))); 

	////special case 3
	//curve_list.push_back(Curve_2(Point_2(3,3), Point_2(3,5))); 
	//curve_list.push_back(Curve_2(Point_2(3,5), Point_2(3,6))); 
	//curve_list.push_back(Curve_2(Point_2(3,6), Point_2(3,7))); 
	//curve_list.push_back(Curve_2(Point_2(3,7), Point_2(11,8))); 
	//curve_list.push_back(Curve_2(Point_2(11,8), Point_2(19,7))); 
	//curve_list.push_back(Curve_2(Point_2(19,7), Point_2(19,6))); 
	//curve_list.push_back(Curve_2(Point_2(19,6), Point_2(19,5))); 
	//curve_list.push_back(Curve_2(Point_2(19,5), Point_2(19,3))); 
	//curve_list.push_back(Curve_2(Point_2(3,3), Point_2(19,3))); 
	//curve_list.push_back(Curve_2(Point_2(3,5), Point_2(19,5))); 
	//curve_list.push_back(Curve_2(Point_2(3,6), Point_2(19,6))); 
	//curve_list.push_back(Curve_2(Point_2(3,7), Point_2(19,7))); 

	////bounding box
	//curve_list.push_back(Curve_2(Point_2(0,0), Point_2(30,0))); 
	//curve_list.push_back(Curve_2(Point_2(0,0), Point_2(0,30))); 
	//curve_list.push_back(Curve_2(Point_2(0,30), Point_2(30,30))); 
	//curve_list.push_back(Curve_2(Point_2(30,30), Point_2(30,0))); 


	// create points

	//house
	//const int num_of_points = 10;
	//Point_2 *points[num_of_points];

	//Point_2 p0(1,2);                       
	//Point_2 p1(2,6); 
	//Point_2 p2(3,2); 
	//Point_2 p3(2,3);                          
	//Point_2 p4(0,2); 
	//Point_2 p5(2,4);
	//Point_2 p6(2,2);                          
	//Point_2 p7(4,4); 
	//Point_2 p8(2,8);
	//Point_2 p9(6,2);

	//points[0] = &p0;
	//points[1] = &p1;
	//points[2] = &p2;
	//points[3] = &p3;
	//points[4] = &p4;
	//points[5] = &p5;
	//points[6] = &p6;
	//points[7] = &p7;
	//points[8] = &p8;
	//points[9] = &p9;

	//int start = 9;
	//int end = num_of_points;

	////specail cases
	//const int num_of_points = 3;
	//Point_2 *points[num_of_points];

	//Point_2 p0(11,2);                       
	//Point_2 p1(11,9); 
	//Point_2 p2(11,4); 

	//points[0] = &p0;
	//points[1] = &p1;
	//points[2] = &p2;

	//int start, end; 
	////specail case 1 
	////start = 0;  end = 1;

	////specail case 2 
	////start = 1;  end = 2;

	////all cases
	//start = 0;	end = num_of_points;
	

	/***************************
	const char * out_filename = "landmarks.out";
	const char * in_filename = "landmarks.in";
	
	std::ofstream outfile(out_filename);
	if (!outfile.is_open()) 
	{
		std::cerr << "Cannot open file " << out_filename << "!" << std::endl;
		return -1;
	}

	std::ifstream infile(in_filename);
	if (!infile.is_open()) 
	{
	std::cerr << "Cannot open file " << in_filename << "!" << std::endl;
	return -1;
	}

	//write and then read
	outfile.clear();
	outfile << pm;
	pm.read(infile);
	***************************/ 
/*
const int num_of_points = 1;
	Point_2 *points[num_of_points];
	int start, end; 
	start = 0;	end = num_of_points;

	//read the point to locate from a file
	std::ifstream inp_pnt_file(inp_pnt_filename);
	if (!inp_pnt_file.is_open()) 
	{
		std::cerr << "Cannot open file " << inp_pnt_filename << "!" << std::endl;
		getchar(); return -1;
	}

	// read
    Coord_type x, y;
    inp_pnt_file >> x >> y;
	Point_2 p0(x,y);     
	points[0] = &p0;

	//open output file
	std::ofstream outfile(out_filename);
	if (!outfile.is_open()) 
	{
		std::cerr << "Cannot open file " << out_filename << "!" << std::endl;
		return -1;
	}

	//locate
	for( int point_index = start; point_index < end; point_index++)
	{
		Point_2 *pntt = points[point_index];
		Point_2 pnt = *pntt;
		std::cout << std::endl << std::endl;
		std::cout << "-------- point number "<< point_index <<" is " << pnt << std::endl << std::endl;

		//locate point
		e = pm.locate(pnt,lt);

		//print output
		if (lt==Pm::UNBOUNDED_FACE) std::cout << "Unbounded face" << std::endl;
		else if (lt==Pm::FACE) std::cout << "Face that is left of :  ";
		else if (lt==Pm::EDGE) std::cout << "EDGE : " ;
		else if (lt==Pm::VERTEX) std::cout << "VERTEX : " ;
		else std::cout << "Unknown locate type" << std::endl;
		//print e
		if (lt != Pm::UNBOUNDED_FACE)
		{
			std::cout << "e = "<< e->source()->point() <<" towards "<< e->target()->point() << std::endl;
			//outfile << e->curve() << std::endl;
			//if (lt == Pm::FACE) {
			//	Ccb_halfedge_circulator curr = e;
			//	curr ++;
			//	while (curr != e) {
			//		outfile << curr->curve() << std::endl;	
			//		curr ++;
			//	}
			//}
		}
	}
*/
