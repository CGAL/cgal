/*! \file
 * A construction of a quadrilateral-shaped Planar_map.
 * In this example a quadrilateral-shaped Planar_map that uses the naive point
 * location strategy is constructed utilising various edge-insertion
 * operations, that is insert_in_face_interior(), insert_from_vertex(), and
 * insert_at_vertices(). The later is verfiried to return the halfedge that is
 * incident to the unbounded face.
 */

#include "short_names.h"

#include <CGAL/Homogeneous.h>
#include <CGAL/Pm_segment_traits_2.h>
#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/Pm_naive_point_location.h>
#include <iostream>

typedef CGAL::Homogeneous<long>           Kernel;
typedef CGAL::Pm_segment_traits_2<Kernel> Traits;
typedef Traits::Point_2                   Point_2;
typedef Traits::X_monotone_curve_2        X_monotone_curve_2;
typedef CGAL::Pm_default_dcel<Traits>     Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>   Planar_map;

int main()
{
  // Create an instance of a Planar_map with a naive point location strategy:
  CGAL::Pm_naive_point_location<Planar_map> naive_pl;
  Planar_map pm(&naive_pl);
  X_monotone_curve_2 cv[4];
  int i;

  CGAL::set_ascii_mode(std::cout);

  Point_2 a1(0, 0), a2(1, 0), a3(1, 1), a4(0, 1);

  // Create the curves:
  cv[0] = X_monotone_curve_2(a1, a2);
  cv[1] = X_monotone_curve_2(a2, a3);
  cv[2] = X_monotone_curve_2(a3, a4);
  cv[3] = X_monotone_curve_2(a4, a1);

  // Insert the curves into the Planar_map:
  Planar_map::Halfedge_handle e[4];  

  e[0] = pm.insert_in_face_interior(cv[0], pm.unbounded_face());
  for (i = 1; i < 3; i++)
    e[i] = pm.insert_from_vertex(cv[i], e[i-1]->target());
  e[3] = pm.insert_at_vertices(cv[3], e[0]->source(), e[2]->target());

  // Check for correct topology:
  CGAL_assertion(e[3]->face()->is_unbounded());

  // Check the validity of the map:
  CGAL_assertion(pm.is_valid());
    
  return 0;  
}
