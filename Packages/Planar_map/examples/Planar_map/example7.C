/*! \file
 * A construction of a Planar_map templated with staight traits.
 * In this example a Planar_map that uses the default point location strategy
 * is constructed. The map is templated with the straight traits which is
 * templated with the Cartesian coordinates over doubles in turn. The example
 * demonstrates insertions and vertical ray shooting.
 */

#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Pm_straight_traits.h>
#include <CGAL/IO/Straight_2_stream.h>
#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/Pm_walk_along_line_point_location.h>
#include <iostream>
#include <iterator>
#include <algorithm>

typedef CGAL::Quotient<CGAL::MP_Float>          Number_type;
typedef CGAL::Cartesian<Number_type>            Kernel;
typedef CGAL::Pm_straight_traits<Kernel>        Traits;
typedef Traits::Point_2                         Point_2;
typedef Traits::X_curve_2                       X_curve_2;
typedef Traits::X_bounded_curve                 Segment;
typedef CGAL::Pm_default_dcel<Traits>           Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>         Planar_map;

CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(X_curve_2);

int main() 
{
  // Create an instance of a planar map:
  CGAL::Pm_walk_along_line_point_location<Planar_map> walk_pl;
  Planar_map pm(&walk_pl);
  X_curve_2 cv[5];

  Point_2 a1(100, 0), a2(20, 50), a3(180, 50), a4(100, 100);

  // Create the curves:
  cv[0] = X_curve_2(Segment(a1, a2));
  cv[1] = X_curve_2(Segment(a1, a3));
  cv[2] = X_curve_2(Segment(a2, a3));
  cv[3] = X_curve_2(Segment(a2, a4));
  cv[4] = X_curve_2(Segment(a3, a4));
  
  std::cout << "The curves of the map :" << std::endl; 
  std::copy(&cv[0], &cv[5], std::ostream_iterator<X_curve_2>(std::cout, "\n"));
  std::cout << std::endl;

  // Insert the five curves into the map:
  std::cout << "Inserting the curves into the map ... ";
  pm.insert(&cv[0], &cv[5]);
  std::cout << ((pm.is_valid()) ? "map valid!" : "map invalid!") << std::endl
            << std::endl;
  
  // Shoot a vertical ray upward from p:
  Point_2 p(95, 30);
  Planar_map::Halfedge_handle e;  
  Planar_map::Locate_type lt;

  std::cout << std::endl << "Upward vertical ray shooting from " << p;
  std::cout << std::endl; 

  e = pm.vertical_ray_shoot(p, lt, true);
  std::cout << "returned the curve " << e->curve() << ", oriented toward " 
            << e->target()->point() << std::endl;
  
  return 0;
}
