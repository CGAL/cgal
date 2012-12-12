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


typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef CGAL::Polyhedron_3<K>                     Polyhedron_3;

typedef K::Segment_3                              Segment_3;
typedef CGAL::Delaunay_triangulation_3<K>         Delaunay;



typedef CGAL::Delaunay_triangulation_sphere_traits_2<K>             Gt;
typedef CGAL::Projection_sphere_traits_3<K>							Gt2;
typedef CGAL::Delaunay_triangulation_sphere_2<Gt>              RTOS;
typedef CGAL::Delaunay_triangulation_sphere_2<Gt2>              RTOS2;
typedef RTOS::Vertex_handle                             Vertex_handle;
typedef RTOS::Face_handle                                 Face_handle;
typedef K::Point_3                                            Point;
typedef RTOS::All_faces_iterator                            Face_iterator;
typedef RTOS::All_vertices_iterator                           Vertex_iterator;
typedef RTOS::Solid_faces_iterator						Solid_faces_iterator;
typedef RTOS::All_edges_iterator						All_edges_iterator;
typedef RTOS::Locate_type                                 Locate_type;
typedef RTOS::Edge                                               Edge;



typedef CGAL::Polyhedron_3<K>                     Polyhedron_3;
typedef K::Segment_3							Segment_3;
typedef CGAL::Delaunay_triangulation_3<K, CGAL::Fast_location> Delaunay_fast;
typedef CGAL::Creator_uniform_3<double, Point>  PointCreator;




int main(){
	
	
	
	int nu_of_pts;
	double radius;
	nu_of_pts =pow(2,10);
	radius=6000000;
	double minDist = radius * pow (2, -25);
	double minDist2 = pow(minDist, 2);
	int invalid = 0;
	CGAL::Timer time;

	CGAL::Random random(nu_of_pts);
	typedef CGAL::Creator_uniform_3<double, Point> Creator;
    CGAL::Random_points_on_sphere_3<Point, Creator> on_sphere(radius);
	
	
	std::vector<Point> points;
	std::vector<Vertex_handle> vertices;
	vertices.reserve(nu_of_pts);
	
	
	for (int count=0; count<nu_of_pts; count++) {
		Point p = *on_sphere;
		points.push_back(p);
		on_sphere++;
	}
	
	
	
	//Triangulation with points on the sphere
	RTOS rtos;
	rtos.set_radius(radius);

	std::cout<<" ***STARTING***"<<std::endl;
	time.start();
	rtos.insert(points.begin(),points.end());
	time.stop();
	std::cout<<"running time triangulation sphere    "<< time.time()<<std::endl;

	

	
	//Triangulation with points on the sphere
	RTOS2 rtos2;
	rtos2.set_radius(radius);
	
	time.reset();
	time.start();
	rtos2.insert(points.begin(),points.end());
	time.stop();
	std::cout<<"running time triangulation sphere projection traits:   "<< time.time()<<std::endl;
	
	
	
		
	Polyhedron_3 poly;
	
	time.reset();
	time.start();
	CGAL::convex_hull_3(points.begin(), points.end(), poly);
	time.stop();
	std::cout << "Static :" << time.time() <<" "<<  std::endl;
	
	poly.clear();
		
	
	
	time.reset();
	time.start();
	CGAL::convex_hull_incremental_3( points.begin(), points.end(), poly, false);
	time.stop();
	std::cout << "incremental EPIC :" << time.time() << std::endl;
	
		
	time.reset();
	time.start();
	Delaunay T_on(points.begin(), points.end());
	time.stop();
	std::cout << "Delaunay on sphere:" << time.time() << std::endl;
	
	time.reset();
	time.start();
	Delaunay_fast T_fast_on(points.begin(), points.end());
	time.stop();
	std::cout << "Delaunay fast location on sphere  :" << time.time() << std::endl;
	
	
	
	
	
}
	
	
	
	
	
	
	
	
	
	
