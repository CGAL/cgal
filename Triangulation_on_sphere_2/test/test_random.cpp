#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_sphere_traits_2.h>
//#include <CGAL/Projection_sphere_traits_3.h>
#include <CGAL/Delaunay_triangulation_sphere_2.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/point_generators_2.h>
#include <fstream>
#include <CGAL/Timer.h>
//#include <CGAL/number_utils.h>
#include <CGAL/squared_distance_3.h>
#include <cmath>

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef CGAL::Delaunay_triangulation_sphere_traits_2<K>             Gt;
//typedef CGAL::Projection_sphere_traits_3<K>							Gt;
typedef CGAL::Delaunay_triangulation_sphere_2<Gt>              RTOS;
typedef RTOS::Vertex_handle                             Vertex_handle;
typedef RTOS::Face_handle                                 Face_handle;
typedef RTOS::Point                                             Point;
typedef RTOS::All_faces_iterator                            Face_iterator;
typedef RTOS::All_vertices_iterator                           Vertex_iterator;
typedef RTOS::Solid_faces_iterator						Solid_faces_iterator;
typedef RTOS::All_edges_iterator						All_edges_iterator;
typedef RTOS::Locate_type                                 Locate_type;
typedef RTOS::Edge                                               Edge;
                              






bool has_face(Face_handle fh, Vertex_handle v0, Vertex_handle v1, Vertex_handle v2){
	bool test1, test2, test3;
	
	
	for(int i=0;i<=2; i++){
		test1 = (v0->point()==fh->vertex(i)->point());
		if(test1) 
			break;
		}
	if(!test1) return false;
	
	for(int i=0;i<=2; i++){
		test2 = v1->point()==fh->vertex(i)->point();
		if(test2) break;
	}
	if(!test2)return false;
	
	for(int i=0; i<=2; i++){
		test3 = v2->point()==fh->vertex(i)->point();
		if(test3)break;
	}	
	if(!test3) return false;
	
	return true;
}
	



bool are_equal(RTOS triA, RTOS triB){
	bool test = false;
	Face_iterator fiA;
	Face_iterator fiB;
	fiA = triA.all_faces_begin();
	fiB = triB.all_faces_begin();
    for( ; fiA != triA.all_faces_end(); ++fiA ){
		//**face of fiA in fiB?
		for( ; fiB != triB.all_faces_end(); ++fiB ){
			test = has_face(fiB, fiA->vertex(0), fiA->vertex(1), fiA->vertex(2));
			if(has_face) break;
			}
		CGAL_assertion(has_face);
		//**	
	}
	return true;
}






int main(){
	int nu_of_pts;
	double radius;
	nu_of_pts =10;
	radius=6000000;
	double minDist = radius * pow (2, -25);
	double minDist2 = pow(minDist, 2);
	int invalid = 0;
	CGAL::Timer t;

	
			
	CGAL::Random random(nu_of_pts);
	typedef CGAL::Creator_uniform_3<double, K::Point_3> Creator;
    CGAL::Random_points_on_sphere_3<K::Point_3, Creator> on_sphere(radius);
	RTOS rtos;
	RTOS rtos2;
	rtos.set_radius(radius);
	rtos2.set_radius(radius);
	
	std::vector<K::Point_3> points;
	std::vector<Vertex_handle> vertices;
		vertices.reserve(nu_of_pts*2);
	
	
	for (int count=0; count<nu_of_pts; count++) {
		K::Point_3 p = *on_sphere;
		points.push_back(p);
		on_sphere++;
	}
	
	
	t.start();
	//====insert new points============
		rtos.insert(points.begin(),points.end());
	//spatial_sort (points.begin(), points.end());
	
	t.stop();
	std::cout<<"running time"<< t.time()<<std::endl;

	std::cout<<"number of vertices    "<<rtos.number_of_vertices()<<std::endl;

	K::Point_3 q = K::Point_3(500,0,0);
	rtos.insert(q);
	rtos.is_valid();
	
	All_edges_iterator eit=rtos.all_edges_begin();
	
	for ( ; eit !=rtos.all_edges_end(); ++eit) 
		CGAL::Object o = rtos.dual(eit);
	
	/*
	
	//*****second triangulation*******
	std::random_shuffle(points.begin(), points.end());
	std::vector<Vertex_handle> vertices2;
	vertices2.reserve(nu_of_pts*2);
	
	for (int count=0; count<nu_of_pts*2; count++) {
		//std::cout<< "================= point number   "<< count+1 <<" =================="<<std::endl;
		K::Point_3 p = points.at(count);
		Vertex_handle v = rtos2.insert(p);
		vertices2.push_back(v);			
		
	}
	rtos2.is_valid();
	
	
			//rtos.show_all();
	t.stop();
		std::cout<<"running time"<< t.time()<<std::endl;
		std::cout<<"number of points"<<std::endl;
	//}
	std::cout<<"comparing"<<std::endl;
	are_equal(rtos, rtos2);
	
	std::cout<<"number of ghost faces  "<<rtos.number_of_ghost_faces()<<std::endl;
	std::cout<<"total faces  "<<rtos.number_of_faces()<<std::endl;
	int count =0;
	
	Solid_faces_iterator sfi= rtos.solid_faces_begin();
	for( ; sfi!= rtos.solid_faces_end(); ++sfi ){
		CGAL_assertion(!sfi->is_ghost());
		count ++;
	}
	std::cout<<"number of solid faces  "<<count<<std::endl;
		
	
	
	
	/*
	 //==remove points=============================
	 std::random_shuffle(vertices.begin(), vertices.end());
	 
	 for( int i=0; i< (int)vertices.size(); i++){
	 rtos.remove(vertices.at(i));
	 std::cout<<rtos.number_of_vertices()<<std::endl;
	 //tos.is_valid();
	 }
	 
	 */
	
	


}