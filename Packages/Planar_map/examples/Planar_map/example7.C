// examples/Planar_map/example7.C
// ------------------------------
#include <CGAL/Cartesian.h>
#include <CGAL/Pm_straight_exact_traits.h>
#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/IO/Straight_2_stream.h>
#include <iostream>

typedef double                              	   Number_type;
typedef CGAL::Cartesian<Number_type>               Coord_t;
typedef CGAL::Pm_straight_exact_traits<Coord_t>    Pmtraits;
typedef Pmtraits::Point                            Point;
typedef Pmtraits::X_curve                          Curve;
typedef Pmtraits::X_bounded_curve                  Segment;
typedef CGAL::Pm_default_dcel<Pmtraits>            Pmdcel;
typedef CGAL::Planar_map_2<Pmdcel,Pmtraits>        Planar_map;


int main() 
{
  // creating an instance of Planar_map_2<pmdcel,pmtraits>
  Planar_map pm;

  Curve cv[5];
  int i;

  Point a1(100, 0), a2(20, 50), a3(180, 50), a4(100, 100);

  // those curves are about to be inserted to pm
  cv[0] = Curve(Segment(a1, a2));
  cv[1] = Curve(Segment(a1, a3));
  cv[2] = Curve(Segment(a2, a3));
  cv[3] = Curve(Segment(a2, a4));
  cv[4] = Curve(Segment(a3, a4));
  
  std::cout << "the curves of the map :" << std::endl; 
  for (i = 0; i < 5; i++)
    std::cout << cv[i] << std::endl;

  std::cout << std::endl;

  // insert the five curves to the map
  std::cout << "inserting the curves to the map..." << std::endl;
  for (i = 0; i < 5; i++)
  {
    std::cout << "inserting curve" << i << std::endl;
    pm.insert(cv[i]);
  }
  
  // check the validity of the map
  std::cout << "check map validity... ";
  if (pm.is_valid())
    std::cout << "map valid!" << std::endl;
  else
    std::cout << "map invalid!" << std::endl;
  std::cout << std::endl;
  
  // vertical ray shooting upward from p
  Point p(95, 30);
  Planar_map::Halfedge_handle e;  
  Planar_map::Locate_type lt;

  std::cout << std::endl << "upward vertical ray shooting from " << p;
  std::cout << std::endl; 

  e=pm.vertical_ray_shoot(p, lt, true);
  std::cout << "returned the curve :" << e->curve() << " oriented toward " 
            << e->target()->point() << std::endl;
  
  return 0;  
}
