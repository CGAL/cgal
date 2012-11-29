#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_sphere_traits_2.h>
#include <CGAL/Triangulation_sphere_2.h>
#include <CGAL/Delaunay_triangulation_sphere_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef CGAL::Delaunay_triangulation_sphere_traits_2<K>             Gt;
typedef CGAL::Delaunay_triangulation_sphere_2<Gt>                 DTOS;
typedef DTOS::Point												Point;

template<class DTOS>
void test (){
	DTOS dtos;
	double radius =1000;
	dtos.set_radius(radius);
	std::vector<Point> points;
//legal points
	Point p1 = Point (radius/sqrt(2), 0, radius/sqrt(2));
	points.push_back(p1);
	
	
	Point p2 = Point (-radius/sqrt(3), radius/sqrt(3),-radius/sqrt(3));
	points.push_back(p2);
	
	
	Point p3 = Point ( radius,0,0);
	points.push_back(p3);
		
	Point p4 = Point (-radius/sqrt(2), -radius/sqrt(2) , 0);
	points.push_back(p4);
	
//points not on sphere
	Point p5 = Point (-1.5*radius,0,0);
	points.push_back(p5);
	
	
	Point p6 = Point (radius, 0, -0.5 * radius);
	points.push_back(p6);
	
	
//points to close
	Point p7 = Point (radius/sqrt(2), radius*pow(2, -15), radius/sqrt(2));
	points.push_back(p7);
		
	Point p8 = Point(radius,0, radius *pow(2, -14));
	points.push_back(p8);
	
	
						 
	dtos.insert(points.begin(), points.end());
	dtos.is_valid();
	assert(dtos.number_of_vertices()==4);
	
    dtos.clear();
	dtos.set_radius(radius);
	dtos.insert(p1);
	assert(dtos.number_of_vertices()==1);
	
	dtos.insert(p1);
	assert(dtos.number_of_vertices()==1);
	
	dtos.insert(p2);
	assert(dtos.number_of_vertices()==2);
	dtos.insert(p2);
	assert(dtos.number_of_vertices()==2);
	
	dtos.insert(p3);
	assert(dtos.number_of_vertices()==3);
	dtos.insert(p3);
	assert(dtos.number_of_vertices()==3);
	
	dtos.insert(p4);
	assert(dtos.number_of_vertices()==4);
	dtos.insert(p4);
	assert(dtos.number_of_vertices()==4);
	
	dtos.insert(p5);
	assert(dtos.number_of_vertices()==4);

	dtos.insert(p6);
	assert(dtos.number_of_vertices()==4);

	dtos.insert(p7);
	assert(dtos.number_of_vertices()==4);

	dtos.insert(p8);
	assert(dtos.number_of_vertices()==4);

	dtos.insert(p8);
	assert(dtos.number_of_vertices()==4);
	
	
}





int main(){
	std::cout<<"Test inserting illegal points"<<std::endl;
	test<DTOS>();
	
	return 0;
}