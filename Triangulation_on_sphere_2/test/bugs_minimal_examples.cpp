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
	typedef K::Point_3 Point_3;
	//DTOS_exact dtos;
	RTOS rtos;
	
	Point_3 p1=Point_3(1/sqrt(2), 1/sqrt(2), 0);
	Point_3 p2 = Point_3(-1/sqrt(2), 1/sqrt(2), 0);
	Point_3 p3 = Point_3(-1/sqrt(2), -1/sqrt(2), 0);
	Point_3 p4 = Point_3(1/sqrt(2), -1/sqrt(2), 0);
	Point_3 p5 = Point_3(1,0,0);
  
		Vertex_handle v1 = rtos.insert(p1);
		Vertex_handle v2 = rtos.insert(p2);
		Vertex_handle v3 = rtos.insert(p3);
		Vertex_handle v4 = rtos.insert(p4); 
	rtos.is_valid();
	     Vertex_handle v5 = rtos.insert(p5);
	rtos.show_all();
	std::cout<<"dimension" << rtos.dimension()<<std::endl;
	
	rtos.remove(v1);
	std::cout<<"dimension" << rtos.dimension()<<std::endl;
	rtos.is_valid();
	
	rtos.remove(v2);
	std::cout<<"dimension" << rtos.dimension()<<std::endl;
	rtos.is_valid();
	
	rtos.remove(v3);
	std::cout<<"dimension" << rtos.dimension()<<std::endl;
	rtos.is_valid();
	
	rtos.remove(v4);
	std::cout<<"dimension" << rtos.dimension()<<std::endl;
	rtos.is_valid();
	
	//rtos.remove(v5);
	std::cout<<"dimension" << rtos.dimension()<<std::endl;
	rtos.is_valid();
	
	
		 
	}

