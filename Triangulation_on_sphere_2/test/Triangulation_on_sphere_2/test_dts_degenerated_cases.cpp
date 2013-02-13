#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_sphere_traits_2.h>
#include <CGAL/Projection_sphere_traits_3.h>
#include <CGAL/Triangulation_sphere_2.h>
#include <CGAL/Delaunay_triangulation_sphere_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef CGAL::Delaunay_triangulation_sphere_traits_2<K>             Gt;
typedef CGAL::Projection_sphere_traits_3<K>						  Projection_traits;
typedef CGAL::Delaunay_triangulation_sphere_2<Gt>                 DTOS;
typedef CGAL::Delaunay_triangulation_sphere_2<Projection_traits>      PDTOS;
//typedef DTOS::Point												Point;
typedef DTOS::Face_handle										Face_handle;



template < class Face_handle, class Vertex_handle>
bool has_face(const Face_handle& fh,const Vertex_handle& v0,const Vertex_handle& v1, const Vertex_handle& v2){
//bool has_face(typename Triangul::Face_handle fh,typename Triangul::Vertex_handle v0, typename Triangul::Vertex_handle v1,typename Triangul::Vertex_handle v2){ 

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

template <class Triangul>
bool are_equal(Triangul triA, Triangul triB){
	typedef typename Triangul::All_faces_iterator                            Face_iterator;
	typedef typename Triangul::All_vertices_iterator                          Vertex_iterator;
	typedef typename Triangul::Face_handle    Face_handle;
	typedef typename Triangul::Vertex_handle     Vertex_handle;
	typedef  typename  Triangul::Point           Point;

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
	 for( ; fiA != triA.all_faces_end(); ++fiA ){
		found=false;
		for(fiB=triB.all_faces_begin(); fiB!=triB.all_faces_end(); fiB++){
			Face_handle fb = Face_handle(fiB);
			Face_handle fa = Face_handle(fiA);
			Vertex_handle v0 = fa->vertex(0);
			Vertex_handle v1 = fa->vertex(1);
			Vertex_handle v2= fa->vertex(2);
			test = has_face(fb, v0, v1, v2);
			if(test){
				found=true;
				break;
			}
				
		}
		assert(found==true);
		
	}
	if(found)
	  return true;
	return false;
}



//template<class DTOS>
void test_Delaunay(){
	
//tests whether it is possible to insert points in degenerated positions and whether the result is uniquely defined after this.	
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
	
	//Delaunay-triangulation_sphere_traits
	dtos.insert(p1);
	dtos.insert(p2);
	dtos.insert(p3);
	dtos.insert(p4);
	dtos.insert(p5);
	dtos.insert(p6);
	dtos.is_valid();
	
		
	std::random_shuffle(points.begin(), points.end());
	for(int i=0; i<6; i++)
		dtos2.insert(points.at(i));
	dtos2.is_valid();
	
	std::random_shuffle(points.begin(), points.end());
	 for(int i=0; i<6; i++)
	  dtos3.insert(points.at(i));
	dtos3.is_valid();
		
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
	for(int i=0; i<6; i++)
		dtos.insert(points.at(i));
	dtos.is_valid();
	
	
	std::random_shuffle(points2.begin(), points2.end());
	for(int i=0; i<6; i++)
		dtos2.insert(points.at(i));
	dtos2.is_valid();
	
	
	
	std::random_shuffle(points2.begin(), points2.end());
	for(int i=0; i<6; i++)
		dtos3.insert(points.at(i));
	dtos3.is_valid();
	
			  
	assert(are_equal(dtos3, dtos2)==true);		
	assert(are_equal(dtos3, dtos)==true);
	assert(are_equal(dtos, dtos2)==true);	
		

}

void test_Projection(){
	Projection_traits traits(K::Point_3(0,0,0));
	PDTOS pdtos(traits);
	
	Projection_traits::Construct_projected_point_3 cst =
    traits.construct_projected_point_3_object();
	
	
	
	//tests whether it is possible to insert points in degenerated positions and whether the result is uniquely defined after this.	
	double radius = 100;
	double radius2 = radius*radius;
	typedef K::Point_3 Point_3;
	//DTOS dtos;
	PDTOS pdtos2(traits);
	PDTOS pdtos3(traits);
	pdtos.set_radius(radius);
	pdtos2.set_radius(radius);
	pdtos3.set_radius(radius);
	std::vector<K::Point_3> points;
	std::vector<K::Point_3> points2;
	
	
	
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
	
	pdtos.insert(
				boost::make_transform_iterator(points.begin(), cst),
				boost::make_transform_iterator(points.end(), cst)
	);
	pdtos.is_valid();
	
	//Delaunay-triangulation_sphere_traits
		
	
	//insert   coplanar Points. Points are coplanar but not coplanar with the center of the sphere
	pdtos.clear();
	pdtos2.clear();
	pdtos3.clear();
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
	
	pdtos2.insert(
				boost::make_transform_iterator(points.begin(), cst),
				boost::make_transform_iterator(points.end(), cst)
	);
	pdtos2.is_valid();
	
	
	
	
}



int main(){
	std::cout<<"testing with Delaunay_triangulation_sphere_traits:  "<<std::endl;
	//test<DTOS>();
	test_Delaunay();
	std::cout<< "testing with Projection_sphere_traits:  "<<std::endl;
	//test<PDTOS>();
	test_Projection();
return 0;
}