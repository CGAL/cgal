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


typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef CGAL::Regular_triangulation_sphere_traits_2<K>             Gt;
typedef CGAL::Regular_triangulation_on_sphere_2<Gt>              RTOS;





typedef RTOS::Vertex_handle Vertex_handle;

int main(int argc, char* argv[]) 
{
	
	
	double radius = 100;
	double radius2 = radius*radius;
	typedef K::Point_3 Point_3;
	//DTOS_exact dtos;
	RTOS rtos;
	// insert and remove 5 coplanar points. Points are also coplanar with the center of the sphere
	/*Point_3 p1=Point_3(radius/sqrt(2), radius/sqrt(2), 0);
	Point_3 p2 = Point_3(-1*radius/sqrt(2), radius/sqrt(2), 0);
	Point_3 p3 = Point_3(-1*radius/sqrt(2), -1*radius/sqrt(2), 0);
	Point_3 p4 = Point_3(radius/sqrt(2), -1*radius/sqrt(2), 0);
	Point_3 p5 = Point_3(radius,0,0);
	
 Vertex_handle v1 = rtos.insert(p1);
Vertex_handle v2 = rtos.insert(p2);
Vertex_handle v3 = rtos.insert(p3);
Vertex_handle v4 = rtos.insert(p4); 
 Vertex_handle v5 = rtos.insert(p5);
	rtos.is_valid();		
	rtos.remove(v1);
	rtos.remove(v2);
	rtos.remove(v3);
	rtos.remove(v4);
	rtos.is_valid();
	rtos.remove(v5);*/
	
		
	//insert   coplanar Points. Points are coplanar but not coplanar with the center of the sphere
	rtos.clear();
	Point_3 p0 = Point_3(0,0,radius);
	Point_3 p21 = Point_3(1/sqrt(2),1/sqrt(2),sqrt(radius2-1));
	Point_3 p22 = Point_3(-1/sqrt(2), -1/sqrt(2), sqrt(radius2-1));
	Point_3 p23 = Point_3(0,1,sqrt(radius2-1));
	Point_3 p24 = Point_3(1,0,sqrt(radius2-1));
	Point_3 p25 = Point_3(-1/sqrt(2), 1/sqrt(2), sqrt(radius2-1));
	Point_3 p26 = Point_3(1/sqrt(2), -1/sqrt(2), sqrt(radius2-1));
		//Vertex_handle v0 = rtos.insert(p0);
	
	Vertex_handle v21 = rtos.insert(p21);
	Vertex_handle v22 = rtos.insert(p22);
	Vertex_handle v23 = rtos.insert(p23);
	rtos.show_vertex(v23);
	
	
	
	Vertex_handle v24 = rtos.insert(p24); 
	rtos.show_vertex(v24);
	rtos.is_valid();
		Vertex_handle v25 = rtos.insert(p25);
	rtos.is_valid();
	rtos.show_all();
	Vertex_handle v26 = rtos.insert(p26); 			
	rtos.is_valid();
	
	
	
	
	
	//insert 5th point not coplanar
	Point_3 p27 = Point_3(0,0,radius);
	Vertex_handle v27 = rtos.insert(p27);
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
	
	
	
	
	
	
	
		 
	}

