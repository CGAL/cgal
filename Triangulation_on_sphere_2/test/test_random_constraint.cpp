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


//convex Hull
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/algorithm.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/convex_hull_3_to_polyhedron_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/convex_hull_incremental_3.h>
#include <CGAL/Timer.h>
#include <vector>


#include <CGAL/Constrained_Delaunay_triangulation_sphere_2.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef CGAL::Polyhedron_3<K>                     Polyhedron_3;

typedef K::Segment_3                              Segment_3;
typedef CGAL::Delaunay_triangulation_3<K>         Delaunay;



typedef CGAL::Delaunay_triangulation_sphere_traits_2<K>             Gt;
typedef CGAL::Projection_sphere_traits_3<K>							Gt2;
typedef CGAL::Constrained_Delaunay_triangulation_sphere_2<Gt>                   CTOS;
typedef CGAL::Constrained_Delaunay_triangulation_sphere_2<Gt2>                  CTOS2;
typedef K::Point_3                                                 Point;
typedef  CTOS::Vertex_handle								Vertex_handle;




int main(){
	
	
	int nu_of_pts;
	double radius;
	nu_of_pts =pow(2,5);
	radius=sqrt(110);
	CGAL::Timer time;

	CGAL::Random random(nu_of_pts);
	typedef CGAL::Creator_uniform_3<double, Point> Creator;
    CGAL::Random_points_on_sphere_3<Point, Creator> on_sphere(radius);
	
	
	std::vector<Point> points;
	//std::vector<Point> points2(points.size()+1);
	std::vector<Point> points2;
	//std::vector<Vertex_handle> vertices;
	//vertices.reserve(nu_of_pts);
	
	points2.push_back(Point(0,0,0));
	
	for (int count=0; count<nu_of_pts; count++) {
		Point p = *on_sphere;
		points.push_back(p);
		points2.push_back(p);
		on_sphere++;
	}
	//Delaunay_traits
	CTOS rtos;
	rtos.set_radius(radius);

	CTOS cdts;
	cdts.set_radius(10);
	
	
	std::cout<<" ***STARTING***"<<std::endl;
	time.start();
	Point p0 =Point(7, -5, 6);
    Point p3 =Point(-6, 5, -7);
	Point p2 =Point(5, 7, 6);
	Point p1 =Point(5, 6, 7);
	
	 rtos.insert_constraint(p0,p1);
	rtos.insert_constraint(p1,p2);
	rtos.insert_constraint(p2,p3);
	rtos.insert_constraint(p3,p0);

	//rtos.insert_constraint(points.at(4), points.at(5));
	//rtos.insert(points.begin(), points.end());
	
	
	
	//rtos.insert(points.begin(),points.end());
	time.stop();
	rtos.is_valid();
	std::cout<<"triangulation sphere    "<< time.time()<<std::endl;
	std::cout<<"number of vertices   "<<rtos.number_of_vertices()<<std::endl;
	
	
	//assert(rtos.number_of_vertices() == nu_of_pts+4);
	
	
/*
	Point p4 = Point(10/sqrt(2),-10/sqrt(2),0);
	Point p5 = Point (10/sqrt(2), 10/sqrt(2),0);
	Point p6 = Point (10/sqrt(2), 0, 10/sqrt(2));
	Point p7 = Point (10/sqrt(2), 0, -10/sqrt(2));
	
	
	
	//cdts.insert_constraint(p6,p7);
	

	
	cdts.insert(Point(0,0,10));
	cdts.insert(Point(0,0,-10));
	cdts.insert(Point(0,10,0));
	cdts.insert(Point(0,-10,0));
	cdts.insert(Point(10,0,0));
	cdts.insert(Point(-10,0,0));
	cdts.insert_constraint(p4,p5);
	std::cout<<"number of faces   "<< cdts.number_of_faces()<<std::endl;
	cdts.is_valid();
	cdts.show_all();
	
	cdts.insert(p7);
*/	
	
}
	

	
	
	
	
	
	
	
