#include <iostream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_polygon_2.h>
#include <CGAL/Random.h>
#include <CGAL/algorithm.h>

#include <CGAL/double.h>

typedef double RT;

#define VERBOSE

typedef CGAL::Exact_predicates_inexact_constructions_kernel		K;
typedef K::Point_2							Point_2;
typedef K::Vector_2							Vector_2;
typedef K::Triangle_2							Triangle_2;
typedef std::vector<Point_2>						Container;
typedef CGAL::Random_points_in_triangle_2<Point_2>			Point_generator;

const double EPS = 1e-30;

template<class InputIterator>
bool inside_or_close_to_triangle(const Triangle_2& tri,InputIterator begin, InputIterator end) {
	while(begin!=end) {
		K::FT dist = squared_distance(tri,*begin);
		if(dist>EPS) {
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
		Point_2 pts[3];
		for(int j = 0; j < 3; ++j) {
			pts[j]=Point_2(rand.get_double(),rand.get_double());
		}
		Triangle_2 tri(pts[0],pts[1],pts[2]);
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
