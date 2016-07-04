#include <boost/iterator/transform_iterator.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_sphere_traits_2.h>
#include <CGAL/Projection_sphere_traits_3.h>
#include <CGAL/Delaunay_triangulation_sphere_2.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/point_generators_2.h>
#include <fstream>
#include <CGAL/Timer.h>
#include <CGAL/squared_distance_3.h>
#include <cmath>





typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef CGAL::Delaunay_triangulation_sphere_traits_2<K>             Gt;
typedef CGAL::Projection_sphere_traits_3<K>							Gt2;
typedef CGAL::Delaunay_triangulation_sphere_2<Gt>                   RTOS;
typedef CGAL::Delaunay_triangulation_sphere_2<Gt2>                  RTOS2;
typedef K::Point_3                                                 Point;

void test1(){
	RTOS rtos;
	
	rtos.set_radius(10);
	Point a = Point(0,0,10);
	
	Point b = Point (-10/sqrt(3), -10/sqrt(3), 10/sqrt(3));
	Point c = Point (10/sqrt(3), 10/sqrt(3), 10/sqrt(3));
	Point d = Point (10/sqrt(3), -10/sqrt(3), 10/sqrt(3));
	Point e = Point (-10/sqrt(3), 10/sqrt(3), 10/sqrt(3));
	
	
	
	RTOS::Vertex_handle v1 = rtos.insert(a);
	RTOS::Vertex_handle v2 = rtos.insert(b);
	RTOS::Vertex_handle v3 = rtos.insert(c);
	RTOS::Vertex_handle v4 = rtos.insert(d);
	RTOS::Vertex_handle v5 = rtos.insert(e);
	assert(rtos.number_of_ghost_faces()==2);
	assert(rtos.dimension()==2);
	
  	rtos.remove(v1);
	assert(rtos.dimension()==1);
	rtos.remove(v2);
	assert(rtos.dimension()==1);
	rtos.remove(v3);
	assert(rtos.dimension()==0);
	rtos.remove(v4);
	assert(rtos.dimension()==-1);
	rtos.remove(v5);
	assert(rtos.dimension()==-2);
	
}	

void test2(){

	RTOS rtos;
	
	rtos.set_radius(10);
	Point a = Point(0,0,10);	
	Point b = Point (-10/sqrt(3), -10/sqrt(3), 10/sqrt(3));
	Point c = Point (10/sqrt(3), 10/sqrt(3), 10/sqrt(3));
	Point d = Point (10/sqrt(3), -10/sqrt(3), 10/sqrt(3));
	Point e = Point (-10/sqrt(3), 10/sqrt(3), 10/sqrt(3));
	Point f = Point (0,0,-10);
	
	
	RTOS::Vertex_handle v1 = rtos.insert(a);
	RTOS::Vertex_handle v2 = rtos.insert(b);
	RTOS::Vertex_handle v3 = rtos.insert(c);
	RTOS::Vertex_handle v4 = rtos.insert(d);
	RTOS::Vertex_handle v5 = rtos.insert(e);
	RTOS::Vertex_handle v6 = rtos.insert(f);
	
	assert(rtos.number_of_ghost_faces()==0);
	rtos.remove(v6);
	assert(rtos.number_of_ghost_faces()==2);
	assert(rtos.dimension()==2);



}



int main(){
	std::cout<<" test remove with  dim_down"<<std::endl;
	test1();
	test2();
			
	return 0;
}
	
	
	
	
	
	
	
	
	
	
