// file: examples/Arrangement_2/example6.C

#include "short_names.h"

#include <CGAL/Homogeneous.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Arr_2_default_dcel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Pm_walk_along_line_point_location.h>

typedef CGAL::MP_Float                                  NT;
typedef CGAL::Homogeneous<NT>                           Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>              Traits;
typedef Traits::Point_2                                 Point;
typedef Traits::Curve_2                                 Curve;
typedef Traits::X_monotone_curve_2                      X_monotone_curve_2;
typedef CGAL::Arr_2_default_dcel<Traits>                Dcel;
typedef CGAL::Arrangement_2<Dcel,Traits>                Arr_2;

int main() 
{
  // Choose a point location strategy
  CGAL::Pm_walk_along_line_point_location<Arr_2::Planar_map> pl;
  Arr_2 arr(&pl);

  // Insertion of the curves
  arr.insert(Curve(Point(0, 0), Point(3, 3)));
  arr.insert(Curve(Point(0, 3), Point(3, 0))); 
  arr.insert(Curve(Point(0, 1), Point(3, 1)));

  CGAL_assertion(arr.number_of_halfedges() == 18);
  
  // Traversal of the curves
  Arr_2::Curve_iterator cit;
  Arr_2::Edge_iterator eit;
  for (cit = arr.curve_node_begin(); cit != arr.curve_node_end(); ++cit) 
  {
    std::cout << std::endl << "Curve level:" << std::endl << cit->curve() 
	      << std::endl;
    std::cout << "Edge level:" << std::endl;
    for (eit = cit->edges_begin(); eit != cit->edges_end(); ++eit) 
    {
      std::cout << eit->x_curve() << std::endl ;
    }
  }

  return 0;
}
