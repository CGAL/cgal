#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_sphere_traits_2.h>
#include <CGAL/Triangulation_on_sphere_2.h>
#include <CGAL/Delaunay_triangulation_sphere_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef CGAL::Delaunay_triangulation_sphere_traits_2<K>             Gt;
typedef CGAL::Delaunay_triangulation_sphere_2<Gt>                 DTOS;
typedef DTOS::Point												Point;
typedef DTOS::Face_handle                                 Face_handle;
typedef DTOS::Vertex_handle									Vertex_handle;
typedef DTOS::All_faces_iterator                            Face_iterator;
typedef DTOS::All_vertices_iterator                           Vertex_iterator;



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
bool are_equal(DTOS triA, DTOS triB){
	bool test = false;
	bool found = false;
	if (triA.number_of_vertices()!= triB.number_of_vertices())
		return false;
	if (triA.number_of_faces()!= triB.number_of_faces())
		return false;
	if(triA.number_of_ghost_faces()!=triB.number_of_ghost_faces())
		return false;
	Face_iterator fiA;
	Face_iterator fiB;
	fiA = triA.all_faces_begin();
	//fiB = triB.all_faces_begin();
    for( ; fiA != triA.all_faces_end(); ++fiA ){
		found=false;
		//**face of fiA in fiB?
		//for( ; fiB != triB.all_faces_end(); ++fiB ){
		for(fiB=triB.all_faces_begin(); fiB!=triB.all_faces_end(); fiB++){
			test = has_face(fiB, fiA->vertex(0), fiA->vertex(1), fiA->vertex(2));
			//if(has_face) break;
			if(test){
				found=true;
				break;
			}
				
		}
		assert(found==true);
		//**	
	}
	if(found)
	  return true;
	return false;
}



template<class DTOS>
void test(){
	double radius = 100;
	double radius2 = radius*radius;
	typedef K::Point_3 Point_3;
	DTOS dtos;
	DTOS dtos2;
	DTOS dtos3;
	
	dtos.set_radius(radius);
	dtos2.set_radius(radius);
	dtos3.set_radius(radius);
	std::vector<K::Point_3> points;
	std::vector<K::Point_3> points2;
	std::vector<K::Point_3> points3;
	std::vector<K::Point_3> points4;
	
	// insert 5 coplanar points. Points are also coplanar with the center of the sphere
	Point_3 p1=Point_3(radius/sqrt(2), radius/sqrt(2), 0);
	Point_3 p2 = Point_3(-1*radius/sqrt(2), radius/sqrt(2), 0);
	Point_3 p3 = Point_3(-1*radius/sqrt(2), -1*radius/sqrt(2), 0);
	Point_3 p4 = Point_3(radius/sqrt(2), -1*radius/sqrt(2), 0);
	Point_3 p5 = Point_3(radius,0,0);
	Point_3 p6 = Point_3(0,0,radius);
	points.push_back(p1);
	points.push_back(p2);
	points.push_back(p3);
	points.push_back(p4);
	points.push_back(p5);
	points.push_back(p6);
	points.resize(6);
	
	
	dtos.insert(p1);
	dtos.insert(p2);
	dtos.insert(p3);
	dtos.insert(p4);
	dtos.insert(p5);
	dtos.insert(p6);
	dtos.show_all();
	
	/*dtos.insert(points.begin(), points.end());	
	dtos.is_valid();
	dtos.show_all();*/
	
	std::cout<<"TRIANGULATION 1"<<std::endl;
	dtos.show_all();
	
	std::random_shuffle(points.begin(), points.end());
	for(int i=0; i<6; i++)
		dtos2.insert(points.at(i));
	
	
	dtos2.is_valid();
	
	std::cout<<"TRIANGULATION 2"<<std::endl;
	dtos2.show_all();
	
	
	
	std::random_shuffle(points.begin(), points.end());
	 for(int i=0; i<6; i++)
	  dtos3.insert(points.at(i));
	
	dtos3.is_valid();
	
	std::cout<<"TRIANGULATION 3"<<std::endl;
	dtos3.show_all();
	
	assert(are_equal(dtos, dtos2)==true);
	assert(are_equal(dtos3, dtos2)==true);
	assert(are_equal(dtos, dtos3)==true);
	
	
		//insert   coplanar Points. Points are coplanar but not coplanar with the center of the sphere
	dtos.clear();
	dtos2.clear();
	dtos3.clear();
	Point_3 p0 = Point_3(0,0,radius);
	points2.push_back(p0);
	Point_3 p21 = Point_3(1/sqrt(2),1/sqrt(2),sqrt(radius2-1));
	points2.push_back(p21);
	
	Point_3 p22 = Point_3(-1/sqrt(2), -1/sqrt(2), sqrt(radius2-1));
	points2.push_back(p22);
	
	Point_3 p23 = Point_3(0,1,sqrt(radius2-1));
	points2.push_back(p23);
	
	Point_3 p24 = Point_3(1,0,sqrt(radius2-1));
	points2.push_back(p24);
	
	Point_3 p25 = Point_3(-1/sqrt(2), 1/sqrt(2), sqrt(radius2-1));
	points2.push_back(p25);
	
	Point_3 p26 = Point_3(1/sqrt(2), -1/sqrt(2), sqrt(radius2-1));
	points2.push_back(p26);
	
	Point_3 p27 = Point_3(radius, 0 ,0);
	points2.push_back(p27);
	
	
	points2.resize(7);
	dtos.insert(points2.begin(), points2.end());
	dtos.is_valid();
	dtos.show_all();
	
	
	
	
	std::random_shuffle(points2.begin(), points2.end());
	dtos2.insert(points2.begin(), points2.end());
	dtos2.is_valid();
	dtos2.show_all();
	
	
	std::random_shuffle(points2.begin(), points2.end());
	dtos3.insert(points2.begin(), points2.end());
	dtos3.is_valid();
	
			  
	//assert(are_equal(dtos, dtos2)==true);
	assert(are_equal(dtos3, dtos2)==true);		
	assert(are_equal(dtos3, dtos)==true);
	bool test2 = are_equal(dtos2, dtos);
	assert(are_equal(dtos, dtos2)==true);
	
	
	
		

}


int main(){

	test<DTOS>();
return 0;
}