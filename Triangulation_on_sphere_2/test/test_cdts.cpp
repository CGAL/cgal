#include <boost/iterator/transform_iterator.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_sphere_traits_2.h>
#include <CGAL/Projection_sphere_traits_3.h>
#include <CGAL/Delaunay_triangulation_sphere_2.h>
#include <CGAL/point_generators_3.h>
#include <fstream>
#include <CGAL/Timer.h>
#include <CGAL/squared_distance_3.h>
#include <cmath>

#include <CGAL/Constrained_Delaunay_triangulation_sphere_2.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;


typedef CGAL::Delaunay_triangulation_sphere_traits_2<K>             Gt;
typedef CGAL::Projection_sphere_traits_3<K>							Gt2;
typedef CGAL::Constrained_Delaunay_triangulation_sphere_2<Gt>                   CTOS;
typedef CGAL::Constrained_Delaunay_triangulation_sphere_2<Gt2>                  CTOS2;
typedef K::Point_3                                                 Point;
typedef  CTOS::Vertex_handle								Vertex_handle;



template <class Triangulation>
int number_of_constraints(Triangulation tri){
	typedef typename Triangulation::All_edges_iterator Edges_iterator;
	int number_constr =0;
	Edges_iterator ec;
	ec = tri.all_edges_begin();
	for(; ec !=tri.all_edges_end(); ec++){
			if(tri.is_constrained(*ec))
			number_constr ++;
	}
	return number_constr;
}


void test_cdts_1 ( std::vector<Point> points){
	

	int nu_of_pts = points.size()+8;
	CTOS cdts;
	cdts.set_radius(10);
	
	Point p0 = Point (10,0,0);
	Point p1 = Point (-10,0,0);
	Point p2 = Point(0,10,0);
	Point p3 = Point (0,-10,0);
	Point p4 = Point (0,0,10);
	Point p5 = Point (0,0,-10);
	
	Point p6 = Point (10/sqrt(3), 10/sqrt(3), 10/sqrt(3));
	Point p7 = Point (10/sqrt(3), -10/sqrt(3), 10/sqrt(3));
	
	Point p8 = Point(10/sqrt(2), 10/sqrt(2),0);
	
	Vertex_handle v0 = cdts.insert(p0);
	Vertex_handle v1 = cdts.insert(p1);
	Vertex_handle v2 = cdts.insert(p2);
	Vertex_handle v3 = cdts.insert(p3);
	Vertex_handle v4 = cdts.insert(p4);
	Vertex_handle v5 = cdts.insert(p5);
	
	
	
	cdts.insert_constraint(v0, v2);
	cdts.insert_constraint(v0,v4);
	cdts.insert_constraint(v2,v4);
	cdts.insert(p6);
	cdts.insert(p7);
	
	cdts.is_valid();
	assert(cdts.number_of_faces()==12);
	assert(number_of_constraints(cdts)==3);
	
	//add additional random points
	cdts.insert(points.begin(), points.end());
	cdts.is_valid();
	assert(cdts.number_of_vertices()==nu_of_pts);
	assert(number_of_constraints(cdts)==3);
	
	
}


void test_cdts_2 ( std::vector<Point> points){
	
	//"star shaped polygon"
	
	int nu_of_pts = points.size()+10;
    CTOS cdts;
	cdts.set_radius(10);
	
	Point p0 = Point(10/sqrt(3),10/sqrt(3),10/sqrt(3));
	Point p1 = Point(10/sqrt(2),0,10/sqrt(2));
	Point p2 = Point(10/sqrt(3),-10/sqrt(3),10/sqrt(3));
	Point p3 = Point(0, -10/sqrt(2), 10/sqrt(2));
	Point p4 = Point(-10/sqrt(3),-10/sqrt(3),10/sqrt(3));
	Point p5 = Point(-10/sqrt(2), 0, 10/sqrt(2));
	Point p6 = Point(-10/sqrt(3),10/sqrt(3), 10/sqrt(3));
	Point p7 = Point(0, 10/sqrt(2),10/sqrt(2));
	Point p8 = Point (0,0,10);
	Point p9 = Point (0,0,-10);
	/*cdts.insert_constraint(p10,p11);
	cdts.insert_constraint(p11,p12);
	cdts.insert_constraint(p12,p13);
	cdts.insert_constraint(p13,p14);
	cdts.insert_constraint(p14,p15);
	cdts.insert_constraint(p15,p16);
	cdts.insert_constraint(p16,p17);
	cdts.insert_constraint(p17,p10);*/
	
	Vertex_handle v0 = cdts.insert(p0);
	Vertex_handle v1 = cdts.insert(p1);
	Vertex_handle v2 = cdts.insert(p2);
	Vertex_handle v3 = cdts.insert(p3);
	Vertex_handle v4 = cdts.insert(p4);
	Vertex_handle v5 = cdts.insert(p5);
	Vertex_handle v6 = cdts.insert(p6);
	Vertex_handle v7 = cdts.insert(p7);
	Vertex_handle v8 = cdts.insert(p8);
	
	cdts.insert_constraint(v0,v1);
	cdts.insert_constraint(v1,v2);
	cdts.insert_constraint(v2,v3);
	cdts.insert_constraint(v3,v4);
	cdts.insert_constraint(v4,v5);
	cdts.insert_constraint(v5,v6);
	cdts.insert_constraint(v6,v7);
	cdts.insert_constraint(v7,v0);
	
	cdts.insert(p8);
	cdts.insert(p9);
	cdts.is_valid();
	
	
	assert(cdts.number_of_faces()==16);
	cdts.insert(points.begin(), points.end());
	cdts.is_valid();
	int vert3 = cdts.number_of_vertices();
	assert(cdts.number_of_vertices()==nu_of_pts);
}


void test_cdts_3(std::vector<Point> points){
	int nu_of_pts = points.size()+6;
    CTOS cdts;
	cdts.set_radius(10);
	Point p0 = Point(0,0,10);
	Point p1 = Point (0, sqrt(19), 9);
	
	Point p2 = Point(0,0,-10);
	Point p3 = Point(10,0,0);
	
	Point p4 = Point(10/sqrt(3),-10/sqrt(3),10/sqrt(3));
	Point p5 = Point(-10/sqrt(3),-10/sqrt(3),10/sqrt(3));
	
	
	Vertex_handle v0=cdts.insert(p0);
	Vertex_handle v1=cdts.insert(p1);
	Vertex_handle v2=cdts.insert(p2);
	Vertex_handle v3=cdts.insert(p3);
	Vertex_handle v4=cdts.insert(p4);
	Vertex_handle v5=cdts.insert(p5);
	cdts.insert_constraint(v0,v1);
	cdts.insert_constraint(v2,v3);
	cdts.insert_constraint(v4,v5);
	
	
	
	
	/*cdts.insert_constraint(p2,p3);
	cdts.show_all();
	cdts.insert_constraint(p4,p5);
	cdts.insert_constraint(p0,p1);*/
		
	assert(number_of_constraints(cdts)==3);
	
	cdts.insert(points.begin(), points.end());
	assert(cdts.number_of_vertices()==nu_of_pts);
	int test2 = number_of_constraints(cdts);
	assert(number_of_constraints(cdts)==3);
	cdts.is_valid();
	
	cdts.clear();
	cdts.set_radius(10);
	cdts.insert(points.begin(), points.end());
	cdts.insert_constraint(p2,p3);
	cdts.insert_constraint(p0,p1);
	cdts.insert_constraint(p4,p5);
	assert(cdts.number_of_vertices()==nu_of_pts);
	assert(number_of_constraints(cdts)==3);
	cdts.is_valid();
}
	
void test_cdts_4 (std::vector<Point> points){	
	CTOS cdts;
	cdts.set_radius(1);
	double a = 1/sqrt(2);
	Point p0 = Point(0,0,1);
	Point p1 = Point( a,a,0);
	Point p2 = Point( a, -a,0);
	Point p3 = Point( -a, -a,0);
	Point p4 = Point (-a,a,0);
	Point p5 = Point(0,0,-1);
	
	Point p6 = Point(1,0,0);
	//Point p7 = Point(-1,0,0);
	Point p8 = Point(1/sqrt(3), 1/sqrt(3), 1/sqrt(3));
	
	Vertex_handle v0 = cdts.insert(p0);
	Vertex_handle v1 = cdts.insert(p1);
	Vertex_handle v2 = cdts.insert(p2);
	cdts.insert(p3);
	/*cdts.insert_constraint(p0, p1);
	cdts.show_all();
	cdts.insert_constraint(p0, p2);
	//cdts.insert_constraint(p5, p1);
	//cdts.insert_constraint(p5, p2);
	cdts.insert_constraint(p1, p2);*/
	cdts.insert_constraint(v0,v1);
	cdts.insert_constraint(v0,v2);
	cdts.insert_constraint(v1,v2);
	int testa1 = cdts.number_of_vertices();
	int testb1 = number_of_constraints(cdts);
	cdts.insert(p3);
	cdts.insert(p4);
	cdts.insert(p5);
	cdts.insert(p6);
	int testb2 = number_of_constraints(cdts);
	//cdts.insert(p7);
	cdts.insert(p8);
	
	
	int testa = cdts.number_of_vertices();
	int testb = number_of_constraints(cdts);

	
	
	cdts.is_valid();
	assert(number_of_constraints(cdts)==5);
	
	cdts.insert(points.begin(), points.end());
	cdts.is_valid();
	assert(number_of_constraints(cdts)==5);
}
	
/*-----------------------------------------------------------------
 --------------------------main------------------------------------
 */


int main(){
	
	
	int nu_of_pts=10000;
	CGAL::Timer time;
	CTOS cdts;
	CTOS cdts2;
	cdts.set_radius(10);
	cdts2.set_radius(10);
		
	//generating additional randomed points
	CGAL::Random random(nu_of_pts);
	typedef CGAL::Creator_uniform_3<double, Point> Creator;
    CGAL::Random_points_on_sphere_3<Point, Creator> on_sphere(10);
	std::vector<Point> points;
	
	
	for (int count=0; count<nu_of_pts; count++) {
		Point p = *on_sphere;
		points.push_back(p);
		on_sphere++;
	}
	
	
	std::cout<<"starting test_cdts_1  with 3 (connected)constrained edges"<<std::endl;
	test_cdts_1(points);
	
	std::cout<<"starting test_cdts_2 with a star-shaped polygon as constrained input"<<std::endl;
	test_cdts_2(points);
	
	std::cout<<"starting test_cdts_3 with not connected constrained edges"<<std::endl;
	test_cdts_3(points);	
	
	std::cout<<"starting test_cdts_4  insert in constrained edge"<<std::endl;
	test_cdts_4(points);
	
	
	
	
	//not connected constraints
	
	
	
	
	
	
	
		
}
	

	
	
	
	
	
	
	
