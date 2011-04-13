// examples/Planar_map/example1.C
// ------------------------------
#include <CGAL/Cartesian.h>
#include <CGAL/Pm_segment_epsilon_traits.h>
#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>

typedef double                                      Number_type;
typedef CGAL::Cartesian<Number_type>                Coord_t;
typedef CGAL::Pm_segment_epsilon_traits<Coord_t>    Pmtraits;
typedef Pmtraits::Point                             Point;
typedef Pmtraits::X_curve                           Curve;
typedef CGAL::Pm_default_dcel<Pmtraits>             Pmdcel;
typedef CGAL::Planar_map_2<Pmdcel,Pmtraits>         Planar_map;
int main()
{
  // creating an instance of Planar_map_2<Pmdcel,Pmtraits>
  // Pm_naive_point_location_strategy<pmdcel,pmtraits> Pl_strategy;  
  // Planar_map_2<Pmdcel,Pmtraits> pm(&Pl_strategy);
  Planar_map pm;

  Curve cv[5];
  int i;

  Point a1(100, 0), a2(20, 50), a3(180, 50), a4(100, 100);

  // those curves are about to be inserted to pm
  cv[0] = Curve(a1, a2);
  cv[1] = Curve(a1, a3);
  cv[2] = Curve(a2, a3);
  cv[3] = Curve(a2, a4);
  cv[4] = Curve(a3, a4);
  
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
  std::cout << "returned the curve :" << e->curve() <<  " oriented toward " 
            << e->target()->point() << std::endl; 
  return 0;
}

