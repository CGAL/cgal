
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include<CGAL/Constrained_Delaunay_triangulation_sphere_2.h>
#include<CGAL/Delaunay_mesher_2.h>
#include<CGAL/Delaunay_mesh_face_base_2.h>
#include<CGAL/Delaunay_mesh_sphere_size_criteria_2.h>
#include<CGAL/Constrained_triangulation_face_base_sphere_2.h>
#include<CGAL/Delaunay_triangulation_sphere_traits_2.h>
#include<CGAL/Delaunay_mesh_sphere_traits_2.h>
#include <CGAL/point_generators_3.h>
#include<iostream>



//Mesh_2

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel          K;
typedef CGAL::Delaunay_triangulation_sphere_traits_2<K>             Gt;
typedef CGAL::Delaunay_mesh_sphere_traits_2<Gt>						Mgt;
typedef CGAL::Constrained_triangulation_face_base_sphere_2<Mgt> Cfb;
typedef CGAL::Triangulation_vertex_base_2<Mgt> Vb;

typedef CGAL::Delaunay_mesh_face_base_2<Mgt, Cfb> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;

typedef CGAL::Constrained_Delaunay_triangulation_sphere_2<Mgt, Tds> CDTS;

typedef CGAL::Delaunay_mesh_sphere_size_criteria_2<CDTS> Criteria;
typedef CDTS::Vertex_handle Vertex_handle;
typedef CDTS::Point Point;


template <class Output_iterator>
void read_points(const char* file_path,Output_iterator out){
	int nb;
	double Long, lat;
	
	std::ifstream input(file_path);
	if (!input){
		std::cerr << "Error while reading " << file_path << std::endl;
		exit(EXIT_FAILURE);
	}
	
	input >> nb;
	
	for (int i=0;i<nb;++i){
#warning tmp : must handle the sphere on which are points (parameter of the traits?)
		input >> Long;
		input >> lat;
		Long=Long/180.*M_PI;
		lat=lat/180.*M_PI;
		*out++=	Point(100*cos(Long) * cos (lat),	100*sin(Long) * cos (lat),100*sin(lat));    
		
		
	}
}


template <class Output_iterator>
void read_points2(const char* file_path,Output_iterator out){
	int nb;
	double Long, lat;
	
	std::ifstream input(file_path);
	if (!input){
		std::cerr << "Error while reading " << file_path << std::endl;
		exit(EXIT_FAILURE);
	}
	
	input >> nb;
	
	for (int i=0;i<nb;++i){
#warning tmp : must handle the sphere on which are points (parameter of the traits?)
		input >> Long;
		input >> lat;
		/*Long=Long/180.*M_PI;
		lat=lat/180.*M_PI;
		*out++=	Point(100*cos(Long) * cos (lat),	100*sin(Long) * cos (lat),100*sin(lat));   */
		*out++=K::Point_2(Long,lat);
		
		
	}
}





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
	CGAL::refine_Delaunay_mesh_2(cdt, Criteria(1.25, 0.5,cdt.geom_traits()),false);
		
	cdt.is_valid();
}
void test2(double radius){	
	
	
	CDTS cdt;
	double r = radius;
	double a = r/sqrt(3);
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
	CGAL::refine_Delaunay_mesh_2(cdt, Criteria(0.125, 0.5,cdt.geom_traits()),false);
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
	cdt.set_radius(r);
	std::vector<Point> lst_pt;
	char* filename = "/Users/cwerner/CGAL-git/Triangulation_on_sphere_2/test/norway.poly";
	read_points(filename,
				std::back_inserter(lst_pt));	
	
	
	Vertex_handle v;
	cdt.set_radius(100);
	std::vector<Vertex_handle> vertices;
	for(int i=0;i< lst_pt.size(); i++){
		Point p = lst_pt.at(i);
		double x = p.x();
		double y= p.y();
		double z = p.z();
		//std::cout<<lst_pt.at(i)<<std::endl;
		v =cdt.insert(lst_pt.at(i));
		vertices.push_back(v);
	}
	int n=lst_pt.size();
	cdt.insert_constraint(vertices.at(n-1), vertices.at(0));
	int number_vert= cdt.number_of_vertices();
	for (int i=1; i<n; i++)
		cdt.insert_constraint(vertices.at(i), vertices.at(i-1));
	std::cout<< "starting mehsing "<<std::endl;
	int ghosts = cdt.number_of_ghost_faces();
	std::cout<< "ghosts :"<< ghosts<<std::endl;
	std::cout<<"number of vertices before meshing:  "<< cdt.number_of_vertices() <<std::endl;
	cdt.is_valid();
	CGAL::refine_Delaunay_mesh_2(cdt, Criteria(0.125,0,cdt.geom_traits()),false);
	cdt.is_valid();
	std::cout<<"number of vertices after meshing:  "<<cdt.number_of_vertices()<<std::endl;
}
	

void test6(){

	

	
		typedef CGAL::Triangulation_vertex_base_2<K> Vb2;
	typedef CGAL::Delaunay_mesh_face_base_2<K> Fb2;
	typedef CGAL::Triangulation_data_structure_2<Vb2, Fb2> Tds;
	typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT2;
	typedef CGAL::Delaunay_mesh_size_criteria_2<CDT2> Criteria;
	
	typedef CDT2::Vertex_handle Vertex_handle;
	typedef CDT2::Point Point;
	
		
		CDT2 cdt;
		std::vector<K::Point_2> lst_pt;
	char* filename = "/Users/cwerner/CGAL-git/Triangulation_on_sphere_2/test/norway.poly";
	read_points2(filename,
				std::back_inserter(lst_pt));	
	
	
	Vertex_handle v;

	std::vector<Vertex_handle> vertices;
	for(int i=0;i< lst_pt.size(); i++){
		Point p = lst_pt.at(i);
				v =cdt.insert(lst_pt.at(i));
		vertices.push_back(v);
	}
	int n=lst_pt.size();
	cdt.insert_constraint(vertices.at(n-1), vertices.at(0));
	int number_vert= cdt.number_of_vertices();
	for (int i=1; i<n; i++)
		cdt.insert_constraint(vertices.at(i), vertices.at(i-1));
	std::cout<< "starting mehsing 2d"<<std::endl;
	
	std::cout<<"number of vertices before meshing:  "<< cdt.number_of_vertices() <<std::endl;
	cdt.is_valid();
	CGAL::refine_Delaunay_mesh_2(cdt, Criteria(0.125,0.5,cdt.geom_traits()),false);
	cdt.is_valid();
	
	std::cout<<"number of vertices after meshing:  "<<cdt.number_of_vertices()<<std::endl;
	
}


int main()
{
	
	
	/*std::cout<<  "Meshing planar polygon with four points "<<std::endl;
	std::cout<< " radius  : 1"<<std::endl;
	//test1(1);
	
	std::cout<< " radius  : 100"<<std::endl;
	//test1(100);*/

	
	std::cout<<  "Meshing easy polygon  "<<std::endl;
	 std::cout<< " radius  : 1"<<std::endl;
	 test2(1);
	std::cout<< " radius  : 10"<<std::endl;
	test2(10);
		
	std::cout<<"Mehsing star-shaped polygon"<< std::endl;
	
	std::cout<< "radius : 1"<<std::endl;
	test3(1);
	test3(10);
	
	//std::cout<< "Meshing two polygons " <<std::endl;
	//test4(1);
	
	std::cout<<"Meshing Norway"<<std::endl;
	test5(1);
	
	std::cout<<"Meshing Norway"<<std::endl;
	test6();
}
	 
