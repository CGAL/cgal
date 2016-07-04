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
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/algorithm.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/convex_hull_3_to_polyhedron_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/convex_hull_incremental_3.h>
#include <CGAL/Timer.h>
#include <vector>


typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef CGAL::Polyhedron_3<K>                     Polyhedron_3;

typedef K::Segment_3                              Segment_3;
typedef CGAL::Delaunay_triangulation_3<K>         Delaunay;



typedef CGAL::Delaunay_triangulation_sphere_traits_2<K>             Gt;
typedef CGAL::Projection_sphere_traits_3<K>							Gt2;
typedef CGAL::Delaunay_triangulation_sphere_2<Gt>                   RTOS;
typedef CGAL::Delaunay_triangulation_sphere_2<Gt2>                  RTOS2;
typedef K::Point_3                                                 Point;




typedef CGAL::Polyhedron_3<K>                     Polyhedron_3;
typedef K::Segment_3							Segment_3;
typedef CGAL::Delaunay_triangulation_3<K, CGAL::Fast_location> Delaunay_fast;
typedef CGAL::Creator_uniform_3<double, Point>  PointCreator;




int main(){
	
	
	int nu_of_pts;
	double radius;
	nu_of_pts =pow(2,10);
	radius=6000000;
	CGAL::Timer time;

	CGAL::Random random(nu_of_pts);
	typedef CGAL::Creator_uniform_3<double, Point> Creator;
    CGAL::Random_points_on_sphere_3<Point, Creator> on_sphere(radius);
	
	
	std::vector<Point> points;
	std::vector<Point> points2;
		
	points2.push_back(Point(0,0,0));
	
	for (int count=0; count<nu_of_pts; count++) {
		Point p = *on_sphere;
		points.push_back(p);
		points2.push_back(p);
		on_sphere++;
	}
	//Delaunay_traits
	RTOS rtos;
	rtos.set_radius(radius);

	std::cout<<" ***STARTING***"<<std::endl;
	time.start();
	rtos.insert(points.begin(),points.end());
	time.stop();
	std::cout<<"triangulation sphere    "<< time.time()<<std::endl;

	

	
	//Triangulation with points on the sphere (projection_traits)
	Gt2 traits(K::Point_3(0,0,0));
	RTOS2 dtos(traits);
	dtos.set_radius(radius);
	Gt2::Construct_projected_point_3 cst =
    traits.construct_projected_point_3_object();
	
		
	time.reset();
	time.start();
	dtos.insert(
				boost::make_transform_iterator(points.begin(), cst),
				boost::make_transform_iterator(points.end(), cst)
				);
	time.stop();
	std::cout<<"triangulation sphere projection traits:   "<< time.time()<<std::endl;
	
	
	
	
	Polyhedron_3 poly;
	
	time.reset();
	time.start();
	CGAL::convex_hull_3(points2.begin(), points2.end(), poly);
	time.stop();
	std::cout << "Static :" << time.time() <<" "<<  std::endl;

	
	poly.clear();
		
	
	
	time.reset();
	time.start();
	CGAL::convex_hull_incremental_3( points2.begin(), points2.end(), poly, false);
	time.stop();
	std::cout << "incremental EPIC :" << time.time() << std::endl;
	
		
	time.reset();
	time.start();
	Delaunay T;
	T.insert(Point(0,0,0));
    T.insert(points.begin(), points.end());
	time.stop();
	std::cout << "Delaunay on sphere:" << time.time() << std::endl;
	
		
	
	time.reset();
	time.start();
	Delaunay_fast T_fast_on2;
	T_fast_on2.insert(Point(0,0,0));
	T_fast_on2.insert(points.begin(), points.end());
	time.stop();
	std::cout << "Delaunay fast location on sphere  :" << time.time() << std::endl;

	
	
}
	
	
	
	
	
	
	
	
	
	
