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

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Delaunay_triangulation_sphere_traits_2< K > FKD;
template class CGAL::Regular_triangulation_on_sphere_2<FKD>;
typedef CGAL::Regular_triangulation_on_sphere_2<FKD>  DTOS_exact;

typedef CGAL::Regular_triangulation_sphere_traits_2< K > FKR;
template class CGAL::Regular_triangulation_on_sphere_2<FKR>;
typedef CGAL::Regular_triangulation_on_sphere_2<FKR>  RTOS_exact;

int main(int argc, char* argv[]) 
{
	typedef K::Point_3 Point_3;
	DTOS_exact dtos;
  {
		dtos.insert(Point_3(-0.201668, 0.941192, -0.271087));
		dtos.insert(Point_3(0.217998, 0.861414, 0.458741));
		dtos.insert(Point_3(-0.600943, 0.398519, -0.692857));
		dtos.insert(Point_3(-0.134936, 0.465834, 0.874523)); 
		std::cout << "NUMBER OF VERTICES: " << dtos.number_of_vertices() << std::endl;	
  }
  RTOS_exact rtos;
  {
		rtos.insert(Point_3(-0.201668, 0.941192, -0.271087));
		rtos.insert(Point_3(0.217998, 0.861414, 0.458741));
		rtos.insert(Point_3(-0.600943, 0.398519, -0.692857));
		rtos.insert(Point_3(-0.134936, 0.465834, 0.874523)); 
		std::cout << "NUMBER OF VERTICES: " << rtos.number_of_vertices() << std::endl;	
  }

  DTOS_exact dtos2;
  dtos2.insert(Point_3(0.201668, 0.941192, 0.271087));
  dtos2.insert(Point_3(0.217998, 0.861414, 0.458741));
  dtos2.insert(Point_3(0.600943, 0.398519, 0.692857));
  dtos2.insert(Point_3(0.134936, 0.465834, 0.874523));
  std::cout << "NUMBER OF VERTICES: " << dtos2.number_of_vertices()
	    <<'\t'<<dtos2.dimension()<< std::endl;	
  dtos2.insert(Point_3(sqrt(1.0/3.0), sqrt(1.0/3.0), sqrt(1.0/3.0)));
  std::cout << "NUMBER OF VERTICES: " << dtos2.number_of_vertices() << std::endl;	
		
	return 0;
}

