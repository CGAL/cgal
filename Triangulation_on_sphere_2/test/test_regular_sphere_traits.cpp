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
} 
