#include <iostream>
#include <CGAL/point_generators_3.h>
#include <CGAL/Random.h>
#include <CGAL/algorithm.h>
#include <CGAL/double.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <cassert>

typedef double RT;

typedef CGAL::Exact_predicates_inexact_constructions_kernel 		K;
typedef K::Point_3 							Point_3;
typedef K::Vector_3 							Vector_3;
typedef K::Triangle_3 							Triangle_3;
typedef K::Plane_3 							Plane_3;
typedef std::vector<Point_3>						Container;
typedef CGAL::Random_points_in_triangle_3<Point_3> 			Point_generator;
	
const double EPS = 1e-20;

template<class InputIterator>
bool inside_or_close_to_triangle(const Triangle_3& tri,InputIterator begin, InputIterator end) {
	Plane_3 plane = tri.supporting_plane();
	while(begin!=end) {
		K::FT dist = squared_distance(plane,*begin);
		if(dist>EPS) { 
			return false;
		}
		Triangle_3 OAB = Triangle_3(*begin,tri[0],tri[1]);
		Triangle_3 OAC = Triangle_3(*begin,tri[0],tri[2]);
		Triangle_3 OBC = Triangle_3(*begin,tri[1],tri[2]);
		K::FT OAB_area = sqrt(OAB.squared_area());
		K::FT OAC_area = sqrt(OAC.squared_area());
		K::FT OBC_area = sqrt(OBC.squared_area());
		K::FT tri_area = sqrt(tri.squared_area());
		if(fabs(OAB_area+OAC_area+OBC_area-tri_area)>1e-15) {
			return false;
		}
		++begin;
	}
	return true;
}

int main() {
	CGAL::Random rand;
	Container point_set;
	const int number_triangles = 10;
	const int number_points = 1000;
	for(int i = 0; i < number_triangles; ++i) {
		Point_3 pts[3];
		for(int j = 0; j < 3; ++j) {
			pts[j]=Point_3(rand.get_double(),rand.get_double(),rand.get_double());
		}
		Triangle_3 tri(pts[0],pts[1],pts[2]);
		Point_generator g1( pts[0], pts[1], pts[2] );
		Point_generator g2( tri );
		Point_generator g3( g1 );

		point_set.clear();
		CGAL::cpp11::copy_n( g1, number_points,
		               std::back_inserter(point_set));
		assert(inside_or_close_to_triangle(tri,point_set.begin(),point_set.end()));

		point_set.clear();
		CGAL::cpp11::copy_n( g2, number_points,
		               std::back_inserter(point_set));
		assert(inside_or_close_to_triangle(tri,point_set.begin(),point_set.end()));

		point_set.clear();
		CGAL::cpp11::copy_n( g3, number_points,
		               std::back_inserter(point_set));
		assert(inside_or_close_to_triangle(tri,point_set.begin(),point_set.end()));
	}

	return 0;
}
