#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_sphere_traits_2.h>
#include <CGAL/enum.h>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef CGAL::Delaunay_triangulation_sphere_traits_2<K>          Gt;
typedef Gt::Point_2                                             Point;
typedef Gt::Orientation_2                                      Orientation_2;
typedef Gt::Power_test_2                                     Power_test_2;




int main()
{
  Point P0=Point(0,0,0);
  Point P1=Point(1,0,0);
  Point P2=Point(0,1,0); 
  Point P3=Point(0,0,1);
  Point P4=Point(1,1,1);
  Point P5=Point(-1,-1,-1);
  Point P6=Point(0.5,0.5,0);
  Point P7=Point(2,0,0);
  Point P8=Point(0.5,0,0);
  Point P9=Point(1,1,0);
  Point P10=Point(0.25,0.25,0);
  Point P11=Point(0.5,0,0.5);
  Point P12=Point(0.5,0,-0.5);
  Point P13=Point(0,0,0.5);
  Point P14=Point(-1,0,0);

  CGAL::Oriented_side result;

  //default constructor with sphere (0,0,0)
  Gt traits=Gt();

  std::cout<<"Test power_test_2"<<std::endl;
  //Power_test_2(p,q,r,s)
  result=traits.power_test_2_object()(P1,P2,P3,P4);  
  assert(result==CGAL::ON_POSITIVE_SIDE);
  result=traits.power_test_2_object()(P1,P2,P3,P5); 
  assert(result==CGAL::ON_NEGATIVE_SIDE);  
  result=traits.power_test_2_object()(P1,P2,P3,P6);  
  assert(result==CGAL::ON_ORIENTED_BOUNDARY);

  //Power_test_2(p,q,r) where p,q and r are coplanar points with _sphere
  result=traits.power_test_2_object()(P1,P2,P9);
  assert(result==CGAL::ON_POSITIVE_SIDE);
  result=traits.power_test_2_object()(P2,P1,P9);
  assert(result==CGAL::ON_POSITIVE_SIDE);
  result=traits.power_test_2_object()(P1,P2,P6);
  assert(result==CGAL::ON_ORIENTED_BOUNDARY);
  result=traits.power_test_2_object()(P1,P2,P10);
  assert(result==CGAL::ON_NEGATIVE_SIDE);
  result=traits.power_test_2_object()(P2,P1,P10);
  assert(result==CGAL::ON_NEGATIVE_SIDE);

  //power_test_2(p,q) where p, q and sphere are colinear
  result=traits.power_test_2_object()(P1,P7);
  assert(result==CGAL::ON_POSITIVE_SIDE);
  result=traits.power_test_2_object()(P1,P1);
  assert(result==CGAL::ON_ORIENTED_BOUNDARY);
  result=traits.power_test_2_object()(P1,P8); 
  assert(result==CGAL::ON_NEGATIVE_SIDE);  


  //with other sphere than (0,0,0)
  Gt traits2=Gt(Point(0,0,1));
  //Power_test_2(p,q,r,s)
  result=traits2.power_test_2_object()(P0,P2,P1,P5);  
  assert(result==CGAL::ON_POSITIVE_SIDE);
  result=traits2.power_test_2_object()(P0,P2,P1,P4); 
  assert(result==CGAL::ON_NEGATIVE_SIDE);  
  result=traits2.power_test_2_object()(P0,P1,P1,P10);  
  assert(result==CGAL::ON_ORIENTED_BOUNDARY);

  //Power_test_2(p,q,r) where p,q and r are coplanar points
  result=traits2.power_test_2_object()(P14,P1,P12);
  assert(result==CGAL::ON_POSITIVE_SIDE);
  result=traits2.power_test_2_object()(P1,P2,P6);
  assert(result==CGAL::ON_ORIENTED_BOUNDARY);
  result=traits2.power_test_2_object()(P14,P1,P11);
  assert(result==CGAL::ON_NEGATIVE_SIDE);

  //power_test_2(p,q) where p, q and sphere are colinear
  result=traits2.power_test_2_object()(P13,P0);
  assert(result==CGAL::ON_POSITIVE_SIDE);
  result=traits2.power_test_2_object()(P13,P13);
  assert(result==CGAL::ON_ORIENTED_BOUNDARY);
  result=traits2.power_test_2_object()(P0,P13); 
  assert(result==CGAL::ON_NEGATIVE_SIDE);
	
	
	std::cout<<"Test Orientation_2"<<std::endl;	
	Point p21 = Point(0.5, 0.5, sqrt(0.75));
	Point p22 = Point(- 0.5, 0.5, sqrt(0.75));
	Point p23 = Point(0, -1, sqrt(0.75));
	Point p24 = Point(0,0,-1);
	Point p25 = Point(0,0,-1.5);
	Point p26 = Point(1,1,2*sqrt(0.75));
	
	result = traits.orientation_2_object()(p21,p22,p23);
	assert (result == CGAL::ON_POSITIVE_SIDE);
	result = traits.orientation_2_object()(p23,p22,p21);
	assert (result == CGAL::ON_NEGATIVE_SIDE);
	result = traits.orientation_2_object()(p24,p21,p22, p23);
	assert (result == CGAL::ON_POSITIVE_SIDE);
	
	std::cout<<"Test Coradial_sphere_2"<<std::endl;
	bool coradial;
	coradial = traits.coradial_sphere_2_object()(p24,p25);
	assert (coradial);
	coradial = traits.coradial_sphere_2_object()(p25, p24);
	assert(coradial);
	coradial = traits.coradial_sphere_2_object()(p21, p26);
	assert(coradial);
	coradial = traits.coradial_sphere_2_object()(p22, p26);
	assert(!coradial);

	std::cout<<"Test Inside_cone_2"<<std::endl;
	Point p31 =  Point(1,1,1);
	Point p32 = Point(-1,1,1);
	Point p33 = Point(0, sqrt(1.5),sqrt(1.5));
	Point p34 = Point(0, -sqrt(2),-1);
	Point p35= Point(0.9,0.9,0.9);
	Point p36 = Point(0, sqrt(2.5),0);
	
	bool inside;
	//conde defined by 0, p31 and p32
	inside = traits.inside_cone_2_object()(p31,p32,p33);
	assert(inside);
	inside = traits.inside_cone_2_object()(p31,p32,p34);
	assert(!inside);
	inside = traits.inside_cone_2_object()(p31,p32,p35);
	assert(!inside);
	//not coplanar
	inside = traits.inside_cone_2_object()(p31,p32,p36);
	assert(!inside);
	
	
	
	
	
	
	
} 
