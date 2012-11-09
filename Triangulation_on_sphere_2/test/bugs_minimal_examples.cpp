//#define CGAL_NO_STATIC_FILTERS

#include <CGAL/basic.h>

#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Gmpq.h>

#include <CGAL/Delaunay_triangulation_sphere_traits_2.h>
#include <CGAL/Delaunay_triangulation_sphere_filtered_traits_2.h>

#include <CGAL/Regular_triangulation_sphere_traits_2.h>
#include <CGAL/Regular_triangulation_on_sphere_2.h>

#include <CGAL/Random.h>
#include <CGAL/Timer.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <iostream>
#include <fstream>
#include <string>
#include <CGAL/Triangulation_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef CGAL::Regular_triangulation_sphere_traits_2<K>             Gt;
typedef CGAL::Regular_triangulation_on_sphere_2<Gt>              RTOS;
//typedef CGAL::Triangulation_2<Gt> RTOS;
typedef RTOS::Vertex_handle                             Vertex_handle;
typedef RTOS::Face_handle                                 Face_handle;
typedef RTOS::Point                                             Point;
typedef RTOS::Faces_iterator                            Face_iterator;
typedef RTOS::Vertices_iterator                           Vertex_iterator;
typedef RTOS::Locate_type                                 Locate_type;
typedef RTOS::Edge                                               Edge;



typedef RTOS::Vertex_handle Vertex_handle;


bool has_face(Face_handle fh, Vertex_handle v0, Vertex_handle v1, Vertex_handle v2){
	bool test1, test2, test3;
	
	
	for(int i=0;i<=2; i++){
		std::cout<<v0->point()<<std::endl;
		std::cout<<fh->vertex(i)->point()<<std::endl;
		test1 = (v0->point()==fh->vertex(i)->point());
		if(test1) 
			break;
	}
	if(!test1) return false;
	
	for(int i=0;i<=2; i++){
		test2 = v1->point()==fh->vertex(i)->point();
		if(test2) break;
	}
	if(!test2)return false;
	
	for(int i=0; i<=2; i++){
		test3 = v2->point()==fh->vertex(i)->point();
		if(test3)break;
	}	
	if(!test3) return false;
	
	return true;
}
bool are_equal(RTOS triA, RTOS triB){
	bool test = false;
	Face_iterator fiA;
	Face_iterator fiB;
	fiA = triA.faces_begin();
	fiB = triB.faces_begin();
    for( ; fiA != triA.faces_end(); ++fiA ){
		//**face of fiA in fiB?
		for( ; fiB != triB.faces_end(); ++fiB ){
			test = has_face(fiB, fiA->vertex(0), fiA->vertex(1), fiA->vertex(2));
			if(has_face) break;
		}
		CGAL_assertion(has_face);
		//**	
	}
	return true;
}














int main(int argc, char* argv[]) 
{
	
	
	double radius = 100;
	double radius2 = radius*radius;
	typedef K::Point_3 Point_3;
	//DTOS_exact dtos;
	RTOS rtos;
	RTOS rtos2;
	RTOS rtos3;
	
	std::vector<K::Point_3> points;
	std::vector<K::Point_3> points2;
/*
	// insert and remove 5 coplanar points. Points are also coplanar with the center of the sphere
		Point_3 p1=Point_3(radius/sqrt(2), radius/sqrt(2), 0);
	Point_3 p2 = Point_3(-1*radius/sqrt(2), radius/sqrt(2), 0);
	Point_3 p3 = Point_3(-1*radius/sqrt(2), -1*radius/sqrt(2), 0);
	Point_3 p4 = Point_3(radius/sqrt(2), -1*radius/sqrt(2), 0);
	Point_3 p5 = Point_3(radius,0,0);
	Point_3 p6 = Point_3(0,0,radius);
	points.push_back(p1);
	points.push_back(p2);
	points.push_back(p3);
	points.push_back(p4);
	points.push_back(p5);
	points.push_back(p6);
	points.resize(5);

 Vertex_handle v1 = rtos.insert(p1);
	Vertex_handle v4 = rtos.insert(p4);
Vertex_handle v2 = rtos.insert(p2);
	Vertex_handle v5 = rtos.insert(p5);
Vertex_handle v3 = rtos.insert(p3);
	Vertex_handle v6 = rtos.insert(p6);

 
	rtos.is_valid();

	
std::random_shuffle(points.begin(), points.end());
	for(int i=0; i<5; i++){
		std::cout<<points.at(i)<<std::endl;
		rtos2.insert(points.at(i));
	}
	
	are_equal(rtos, rtos2);
	//rtos.remove(v1);
	//rtos.remove(v2);
	//rtos.remove(v3);
	//rtos.remove(v4);
	//rtos.is_valid();
	//rtos.remove(v5);
	*/
	//insert   coplanar Points. Points are coplanar but not coplanar with the center of the sphere
	rtos.clear();
	Point_3 p27 = Point_3(0,0,radius);
	points2.push_back(p27);
	Point_3 p0 = Point_3(0,0,radius);
	points2.push_back(p0);
	Point_3 p21 = Point_3(1/sqrt(2),1/sqrt(2),sqrt(radius2-1));
	points2.push_back(p21);
	
	Point_3 p22 = Point_3(-1/sqrt(2), -1/sqrt(2), sqrt(radius2-1));
	points2.push_back(p22);
	
	Point_3 p23 = Point_3(0,1,sqrt(radius2-1));
	points2.push_back(p23);
	
	Point_3 p24 = Point_3(1,0,sqrt(radius2-1));
	points2.push_back(p24);
	
	Point_3 p25 = Point_3(-1/sqrt(2), 1/sqrt(2), sqrt(radius2-1));
	points2.push_back(p25);
	
	Point_3 p26 = Point_3(1/sqrt(2), -1/sqrt(2), sqrt(radius2-1));
	points2.push_back(p26);
		//Vertex_handle v0 = rtos.insert(p0);
	
	Vertex_handle v21 = rtos.insert(p21);
	Vertex_handle v22 = rtos.insert(p22);
	Vertex_handle v23 = rtos.insert(p23);
	rtos.show_vertex(v23);
	
	
	
	Vertex_handle v24 = rtos.insert(p24); 
	Vertex_handle v25 = rtos.insert(p25);
	Vertex_handle v26 = rtos.insert(p26); 			
	Vertex_handle v27 = rtos.insert(p27);
	rtos.is_valid();

	
	
	
	points2.resize(8);
	std::random_shuffle(points2.begin(), points2.end());
	for(int i=0; i<8; i++){
		std::cout<<points2.at(i)<<std::endl;
		rtos2.insert(points2.at(i));
	}
	
	points2.resize(8);
	std::random_shuffle(points2.begin(), points2.end());
	for(int i=0; i<8; i++){
		std::cout<<points2.at(i)<<std::endl;
		rtos3.insert(points2.at(i));
	}
	
	
	
	are_equal(rtos3,rtos2);
	/*	
	
	
	//insert 5th point not coplanar
	//Point_3 p27 = Point_3(0,0,radius);
	//Vertex_handle v27 = rtos.insert(p27);
	rtos.show_all();
	rtos.is_valid();
	
	
	rtos.remove(v21);
	rtos.is_valid();
	rtos.show_all();
	
	rtos.remove(v22);
	rtos.is_valid();
	rtos.show_all();
	
	rtos.remove(v23);
	rtos.is_valid();
	rtos.show_all();
	
	
	//3 edges and faces left
	rtos.show_vertex(v24);	
	rtos.remove(v24);
	rtos.is_valid();
	rtos.show_all();
	
	rtos.remove(v25);
	rtos.is_valid();
	rtos.show_all();
	
	rtos.remove(v26);
	rtos.is_valid();
	rtos.show_all();
	
	
	rtos.remove(v27);
	
	*/
	
	
	
	
	
		 
	}

