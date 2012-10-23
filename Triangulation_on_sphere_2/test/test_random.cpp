#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_sphere_traits_2.h>
#include <CGAL/Regular_triangulation_on_sphere_2.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/point_generators_2.h>
#include <fstream>
#include <CGAL/Timer.h>
#include <CGAL/number_utils.h>
#include <CGAL/squared_distance_3.h>
#include <cmath>




typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef CGAL::Regular_triangulation_sphere_traits_2<K>             Gt;
typedef CGAL::Regular_triangulation_on_sphere_2<Gt>              RTOS;
typedef RTOS::Vertex_handle                             Vertex_handle;
typedef RTOS::Face_handle                                 Face_handle;
typedef RTOS::Point                                             Point;
typedef RTOS::Faces_iterator                            Face_iterator;
typedef RTOS::Vertices_iterator                           Vertex_iterator;
typedef RTOS::Locate_type                                 Locate_type;
typedef RTOS::Edge                                               Edge;





bool is_ok(K::Point_3 p, std::vector<K::Point_3> po, double minDist2, int ind)
{
	bool ok= true;
	for(int j= 0; j<ind; j++){
		
		if( squared_distance(p, po.at(j))<=minDist2 ){
		ok = false;
		}
		
	}
	return ok;
}


int main(){
	int nu_of_pts;
	double radius;
	nu_of_pts =10000;
	radius=6000000;
	//radius = 1;
	//double minDist = radius*2e-25;
	double minDist = radius * pow (2, -25);
	double minDist2 = pow(minDist, 2);
	int invalid = 0;
	int random = 0;
	CGAL::Timer t;
	
	
	for(int i=1; i<= 1; i++){
		
		random++;
		std::cout<<" *************************  run  "<< random << "**********************"<<std::endl;
		
	CGAL::Random random(nu_of_pts);
	typedef CGAL::Creator_uniform_3<double, K::Point_3> Creator;
    CGAL::Random_points_on_sphere_3<K::Point_3, Creator> on_sphere(radius);
	RTOS tos;
		
	std::vector<K::Point_3> points;
	std::vector<Vertex_handle> vertices;
		vertices.reserve(nu_of_pts*2);
	
	
	
	for (int count=0; count<nu_of_pts; count++) {
		K::Point_3 p = *on_sphere;
		double tmp = p.x();
		double tmp1 = p.y();
		double tmp2 = p.z();
		
		
		p=K::Point_3(fabs(tmp), fabs(tmp1),  fabs(tmp2));
		//points.push_back(*on_sphere);
		points.push_back(p);
		on_sphere++;
	}
	
	
	for (int count=0; count<nu_of_pts; count++) {
		K::Point_3 p = *on_sphere;
		double tmp = p.x();
		double tmp1 = p.y();
		double tmp2 = p.z();
		
		
		p=K::Point_3(-1*fabs(tmp), -1*fabs(tmp1),  -1*fabs(tmp2));
		
			
		points.push_back(p);
		on_sphere++;
	}
	

		
	t.start();
	//====insert new points============
		
		
	for (int count=0; count<nu_of_pts*2; count++) {
		std::cout<< "================= point number   "<< count+1 <<" =================="<<std::endl;
		K::Point_3 p = points.at(count);
		if(is_ok(p, points, minDist2, count)){
			Vertex_handle v = tos.insert(p);
			vertices.push_back(v);			
		}
		else {
			std::cout<<"invalide point"<<std::endl;
			invalid ++;
		}
		
		//tos.is_valid(true);
	}
	tos.is_valid();
		std::cout<<"starting to remove"<<std::endl;
		
	//==remove points=============================
		std::random_shuffle(vertices.begin(), vertices.end());
	
		for( int i=0; i< (int)vertices.size(); i++){
			tos.remove(vertices.at(i));
			std::cout<<tos.number_of_vertices()<<std::endl;
			//tos.is_valid();
		}
		
		t.stop();
		tos.show_all();
		std::cout<<"running time"<< t.time()<<std::endl;
	}






	
}