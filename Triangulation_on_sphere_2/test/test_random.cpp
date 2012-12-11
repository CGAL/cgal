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
#include <CGAL/algorithm.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>
#include <vector>






typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef CGAL::Delaunay_triangulation_sphere_traits_2<K>             Gt;
typedef CGAL::Projection_sphere_traits_3<K>							Gt2;
typedef CGAL::Delaunay_triangulation_sphere_2<Gt>              RTOS;
typedef CGAL::Delaunay_triangulation_sphere_2<Gt2>              RTOS2;
typedef RTOS::Vertex_handle                             Vertex_handle;
typedef RTOS::Face_handle                                 Face_handle;
typedef RTOS::Point                                             Point;
typedef RTOS::All_faces_iterator                            Face_iterator;
typedef RTOS::All_vertices_iterator                           Vertex_iterator;
typedef RTOS::Solid_faces_iterator						Solid_faces_iterator;
typedef RTOS::All_edges_iterator						All_edges_iterator;
typedef RTOS::Locate_type                                 Locate_type;
typedef RTOS::Edge                                               Edge;
typedef RTOS::Point											Point;


typedef CGAL::Polyhedron_3<K>                     Polyhedron_3;
typedef K::Segment_3							Segment_3;

struct Plane_from_facet {
	Polyhedron_3::Plane_3 operator()(Polyhedron_3::Facet& f) {
		Polyhedron_3::Halfedge_handle h = f.halfedge();
		return Polyhedron_3::Plane_3( h->vertex()->point(),
									 h->next()->vertex()->point(),
									 h->opposite()->vertex()->point());
	}
};



int main(){
	
	
	
	int nu_of_pts;
	double radius;
	nu_of_pts =10;
	radius=6000000;
	double minDist = radius * pow (2, -25);
	double minDist2 = pow(minDist, 2);
	int invalid = 0;
	CGAL::Timer t1, t2, t3;

	CGAL::Random random(nu_of_pts);
	typedef CGAL::Creator_uniform_3<double, Point> Creator;
    CGAL::Random_points_on_sphere_3<Point, Creator> on_sphere(radius);
	
	
	std::vector<Point> points;
	std::vector<Vertex_handle> vertices;
	vertices.reserve(nu_of_pts*2);
	
	
	for (int count=0; count<nu_of_pts; count++) {
		Point p = *on_sphere;
		points.push_back(p);
		on_sphere++;
	}
	//Triangulation with points on the sphere
	RTOS rtos;
	rtos.set_radius(radius);

	std::cout<<" ***STARTING***"<<std::endl;
	t1.start();
	rtos.insert(points.begin(),points.end());
	t1.stop();
	std::cout<<"running time triangulation sphere    "<< t1.time()<<std::endl;

	std::cout<<"number of vertices    "<<rtos.number_of_vertices()<<std::endl;

	
	//Triangulation with points on the sphere
	RTOS2 rtos2;
	rtos2.set_radius(radius);
	
	
	t2.start();
	rtos2.insert(points.begin(),points.end());
	t2.stop();
	std::cout<<"running time triangulation sphere projection traits:   "<< t2.time()<<std::endl;
	
	std::cout<<"number of vertices    "<<rtos2.number_of_vertices()<<std::endl;
	
	rtos.is_valid();
	rtos2.is_valid();
	/*
	//*****Convex Hull***********
		
	// define polyhedron to hold convex hull
	Polyhedron_3 poly;
	
	// compute convex hull of non-collinear points
	t3.start();
	CGAL::convex_hull_3(points.begin(), points.end(), poly);
	t3.stop();
	std::cout<<"running time convex hull t2  "<< t3.time()<<std::endl;
	
	std::cout << "The convex hull contains    " << poly.size_of_vertices() << " vertices" << std::endl;
	
	//*/

}