#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_sphere_traits_2.h>
#include <CGAL/Projection_sphere_traits_3.h>
#include <CGAL/enum.h>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef CGAL::Projection_sphere_traits_3<K>                      Gt;
typedef  Gt::Point_2                                    Point;
typedef  Gt::Orientation_2                              Orientation_2;
typedef  Gt::Power_test_2                               Power_test_2;
typedef  Gt:: Coradial_sphere_2                         Coradial_sphere_2;
typedef  Gt:: Inside_cone_2                             Inside_cone_2;
typedef  Gt::Compute_squared_distance_2                 Compute_squared_distance_2;
typedef  Gt:: Compare_xyz_3                             Compare_xyz_3;


int main()
{
  Gt traits=Gt();
// Operations	
	Power_test_2 power_test = traits.power_test_2_object();
	Orientation_2 orientation = traits.orientation_2_object();
	Coradial_sphere_2 coradial = traits.coradial_sphere_2_object();
	Inside_cone_2 inside_cone = traits.inside_cone_2_object();
	Compute_squared_distance_2 squared_distance= traits.compute_squared_distance_3_object();
	Compare_xyz_3 compare_xyz = traits.compare_xyz_3_object();
	
	traits.set_radius(1);
	
	
	//Testing with Points projected on unit sphere
	
	Point p11 = Point(1.5,1.5,1.5);
	Point p12 = Point(-1,1,1);
	Point p13 = Point(0,0.7,0);
	Point p14 = Point(0,-2,0);
	
	Point p21 = Point(2,0,-1);
	Point p22 = Point(-1.8,0,-0.7);
	Point p23 = Point(0,0,2.6);
	
	
	//Points with same coordinates
	Point p31 = Point(0.6, 0.3,-1);
	Point p32 = Point(0.6, 0.3,-1);
	Point p33 = Point(0.6, 0.3, 1);
	
	//Coradial points
	Point p41 = Point(1,1,1);
	Point p42 = Point(0.7, 0.7, 0.7);
	Point p43 =Point(0.3, 0.6, -2.4);
	Point p44 = Point(0.6, 1.2, -4.8);
	
	
	//Inside Cone
	Point p51 = Point(1,1,1);
	Point p52 = Point (-1, 1, 1);
	Point p53 = Point (1,3,3);//inside
	Point p54 = Point (0.8, 0.8, 0.8); //boundary
	Point p55 = Point(0,-2, -2); //outside
	
	//distance
	Point p61 = Point(1,1,0);
	Point p62 = Point(-1,1,0);
	
	
	std::cout<<"Test Orientation"<<std::endl;
	assert(orientation(p11,p12,p13)==CGAL::NEGATIVE);
	assert(orientation(p13,p12,p11)==CGAL::POSITIVE);
	assert(orientation(p21,p22,p23)==CGAL::ON_ORIENTED_BOUNDARY);
	
	assert(orientation(p11, p12, p13, p14)==CGAL::POSITIVE);
	
	
	
	
	std::cout<<"Test Power_Test_2"<<std::endl;
	assert(power_test(p11, p12, p13, p14)==CGAL::POSITIVE);
	assert(power_test(p14, p11, p12, p13)==CGAL::NEGATIVE);
	assert(power_test(p21,p22,p23,p11)==CGAL::POSITIVE);
	
	
	std::cout<<"Test Coradial_sphere_2"<<std::endl;
	assert(coradial(p41,p41));
	assert(coradial(p41,p42));
	assert(coradial(p43, p44));			  
	
	
	std::cout<<"Test Inside_Cone_sphere_2"<<std::endl;
	assert(inside_cone(p51, p52,p53));
	assert(!inside_cone(p51,p52, p54));	  
	assert(!inside_cone(p51,p52,p55));
		   
	
	std::cout<<"Test Squared_Distance_sphere_2"<<std::endl;
	
	
	
	double dist1 = squared_distance(p41, p42);
	assert(dist1 == 0);
	//double dist2 = squared_distance(p61, p62);
	//assert( dist2 == 2);
	
	
	std::cout<<"Test compare_xyz_3"<<std::endl;
	
	assert(compare_xyz(p31, p31)==CGAL::EQUAL);
	assert(compare_xyz(p31, p32)==CGAL::EQUAL);
	assert(compare_xyz(p31, p33)==CGAL::SMALLER);
	assert(compare_xyz(p33, p31)==CGAL::LARGER);
	
	
	
} 

