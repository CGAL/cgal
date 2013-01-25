
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include<CGAL/Constrained_Delaunay_triangulation_sphere_2.h>
#include<CGAL/Delaunay_mesher_2.h>
#include<CGAL/Delaunay_mesh_face_base_2.h>
#include<CGAL/Delaunay_mesh_size_criteria_2.h>
#include<CGAL/Constrained_triangulation_face_base_sphere_2.h>
#include<CGAL/Delaunay_triangulation_sphere_traits_2.h>
#include<CGAL/Delaunay_mesh_sphere_traits_2.h>
#include <CGAL/point_generators_3.h>
#include<iostream>



typedef CGAL::Exact_predicates_inexact_constructions_kernel          K;
typedef CGAL::Delaunay_triangulation_sphere_traits_2<K>             Gt;
typedef CGAL::Delaunay_mesh_sphere_traits_2<Gt>						Mgt;
typedef CGAL::Constrained_triangulation_face_base_sphere_2<Mgt> Cfb;
typedef CGAL::Triangulation_vertex_base_2<Mgt> Vb;

typedef CGAL::Delaunay_mesh_face_base_2<Mgt, Cfb> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;

typedef CGAL::Constrained_Delaunay_triangulation_sphere_2<Mgt, Tds> CDTS;

typedef CGAL::Delaunay_mesh_size_criteria_2<CDTS> Criteria;
typedef CDTS::Vertex_handle Vertex_handle;
typedef CDTS::Point Point;

void test1(double radius){	
	
		
	CDTS cdt;
	double r = radius;
	double a = r/sqrt(3);
	double b = r/sqrt(2);
	cdt.set_radius(r);
	
	
	
	Vertex_handle va = cdt.insert(Point(a,a,a));
	Vertex_handle vb = cdt.insert(Point(a,-a,a));
	Vertex_handle vc= cdt.insert(Point(-a,-a,a));
	Vertex_handle vd = cdt.insert(Point(-a,a,a));
	
	
	
	cdt.insert_constraint(va, vb);
	cdt.insert_constraint(vb, vc);
	cdt.insert_constraint(vc, vd);
	cdt.insert_constraint(vd, va);
	
	cdt.is_valid();
	CGAL::refine_Delaunay_mesh_2(cdt, Criteria(0.125, 0.5),false);
		
	cdt.is_valid();
}
void test2(double radius){	
	
	
	CDTS cdt;
	double r = radius;
	double a = r/sqrt(3);
	//double a = 0.57735026918962584;
	//double b = 0.707106781187;
	double b = r/sqrt(2);
	cdt.set_radius(r);
	
	
	
	Vertex_handle va = cdt.insert(Point(a,a,a));
	Vertex_handle vb = cdt.insert(Point(a,-a,a));
	Vertex_handle vc= cdt.insert(Point(-r/sqrt(2),-r/sqrt(2),0));
	Vertex_handle vd = cdt.insert(Point(-a,a,a));
	
	
	
	cdt.insert_constraint(va, vb);
	cdt.insert_constraint(vb, vc);
	cdt.insert_constraint(vc, vd);
	cdt.insert_constraint(vd, va);
	
	cdt.is_valid();
	CGAL::refine_Delaunay_mesh_2(cdt, Criteria(0.125, 0.5),false);
	std::cout<<"number of vertices:  "<<cdt.number_of_vertices()<<std::endl;
	cdt.is_valid();
}
void test3(double radius){	
	
	
	CDTS cdt;
	double r = radius;
	double a = r/sqrt(3);
	double b = r/sqrt(2);
	cdt.set_radius(r);
	
	
	
	Vertex_handle va = cdt.insert(Point(b,b,0));
	Vertex_handle vb = cdt.insert(Point(b,0,b));
	Vertex_handle vc = cdt.insert(Point(b,-b,0));
	Vertex_handle vd = cdt.insert(Point(0,-b,b));
	Vertex_handle ve = cdt.insert(Point(-b,-b,0));
	Vertex_handle vf = cdt.insert(Point(-b,0,b));
	Vertex_handle vg= cdt.insert(Point(-b,b,0));
	Vertex_handle vh = cdt.insert(Point(0,b,b));
	
	
	
	cdt.insert_constraint(va, vb);
	cdt.insert_constraint(vb, vc);
	cdt.insert_constraint(vc, vd);
	cdt.insert_constraint(vd, ve);
	cdt.insert_constraint(ve, vf);
	cdt.insert_constraint(vf, vg);
	cdt.insert_constraint(vg, vh);
	cdt.insert_constraint(vh, va);
	
	cdt.is_valid();
	CGAL::refine_Delaunay_mesh_2(cdt, Criteria(0.125, 0.5),false);
	std::cout<<"number of vertices:  "<<cdt.number_of_vertices()<<std::endl;
	cdt.is_valid();
}
void test4(double radius){	
	
	
	CDTS cdt;
	double r = radius;
	double a = r/sqrt(3);
	double b = r/sqrt(2);
	cdt.set_radius(r);
	
	
	
	Vertex_handle va = cdt.insert(Point(b,b,0));
	Vertex_handle vb = cdt.insert(Point(b,-b,0));
	Vertex_handle vc =cdt.insert(Point(0,b,b));
	
	
	Vertex_handle vd = cdt.insert(Point(-b,-b,0));
	Vertex_handle ve= cdt.insert(Point(-b,b,0));
	Vertex_handle vf = cdt.insert(Point(0,-b,b));
	
	
	
	
	cdt.insert_constraint(va, vb);
	cdt.insert_constraint(vb, vc);
	cdt.insert_constraint(va, vc);
	
	cdt.insert_constraint(vd, ve);
	cdt.insert_constraint(ve, vf);
	cdt.insert_constraint(vd, vf);
	
	
	cdt.is_valid();
	CGAL::refine_Delaunay_mesh_2(cdt, Criteria(0.125, 0.5),false);
	std::cout<<"number of vertices:  "<<cdt.number_of_vertices()<<std::endl;
	cdt.is_valid();
}
void test5(double radius){	
	
	
	CDTS cdt;
	double r = radius;
	double a = r/sqrt(3);
	double b = r/sqrt(2);
	cdt.set_radius(r);
	
	
	
	Vertex_handle va = cdt.insert(Point(b,b,0));
	Vertex_handle vb = cdt.insert(Point(b,-b,0));
	Vertex_handle vc =cdt.insert(Point(0,b,b));
	
	
	Vertex_handle vd = cdt.insert(Point(-b,-b,0));
	Vertex_handle ve= cdt.insert(Point(-b,b,0));
	Vertex_handle vf = cdt.insert(Point(0,-b,b));
	
	
	
	
	cdt.insert_constraint(va, vb);
	cdt.insert_constraint(vb, vc);
	cdt.insert_constraint(va, vc);
	
	cdt.insert_constraint(vd, ve);
	cdt.insert_constraint(ve, vf);
	cdt.insert_constraint(vd, vf);
	
	
	cdt.is_valid();
	CGAL::refine_Delaunay_mesh_2(cdt, Criteria(0.125, 0.5),false);
	std::cout<<"number of vertices:  "<<cdt.number_of_vertices()<<std::endl;
	cdt.is_valid();
}
	
	

int main()
{
	/*std::cout<<  "Meshing planar polygon with four points "<<std::endl;
	std::cout<< " radius  : 1"<<std::endl;
	test1(1);
	
	std::cout<< " radius  : 100"<<std::endl;
	test1(100);*/

	
	std::cout<<  "Meshing easy polygon  "<<std::endl;
	 std::cout<< " radius  : 1"<<std::endl;
	 test2(1);
	 
		
	std::cout<<"Mehsing star-shaped polygon"<< std::endl;
	
	std::cout<< "radius : 1<<std::endl"<<std::endl;
	test3(1);
	
	std::cout<< "Meshing two polygons " <<std::endl;
	//test4(1);
}
