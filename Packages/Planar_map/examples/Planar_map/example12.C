/*! \file
 * A construction of a simple planar map.
 * In this example a planar map similar to the planar map of example 1 and 2
 * is constructed.
 */

#include "short_names.h"

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Pm_segment_traits_2.h>
#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/IO/Pm_file_writer.h>
#include <CGAL/IO/write_pm.h>
#include <iostream>

typedef CGAL::Simple_cartesian<float>           Cartesian_kernel;
typedef CGAL::Filtered_kernel<Cartesian_kernel> Kernel;
typedef CGAL::Pm_segment_traits_2<Kernel>       Traits;
typedef Traits::Point_2                         Point_2;
typedef Traits::X_monotone_curve_2              X_monotone_curve_2;
typedef CGAL::Pm_default_dcel<Traits>           Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>         Planar_map;

int main()
{
  // Create an instance of a Planar_map:
  Planar_map pm;

  Point_2 p0(0, 0);
  if (!pm.is_point_in_face(p0, pm.unbounded_face())) {
    std::cerr << "Error: failed to determine point0 in unbounded face!"
              << std::endl;
    return 1;
  }
  
  Point_2 a0(0, 0), a1(2, 0), a2(1, 2);

  // Create the curves:
  X_monotone_curve_2 cv[3];
  cv[0] = X_monotone_curve_2(a0, a1);
  cv[1] = X_monotone_curve_2(a1, a2);
  cv[2] = X_monotone_curve_2(a2, a0);

  std::cout << "The curves of the map :" << std::endl;
  std::copy(&cv[0], &cv[3],
            std::ostream_iterator<X_monotone_curve_2>(std::cout, "\n"));
  std::cout << std::endl;

  // Insert the curves into the Planar_map:
  std::cout << "Inserting the curves to the map ... ";

  Planar_map::Halfedge_handle e[3];  
  e[0] = pm.insert_in_face_interior(cv[0], pm.unbounded_face());
  e[1] = pm.insert_from_vertex(cv[1], e[0]);
  e[2] = pm.insert_at_vertices(cv[2], e[1], e[0]->twin());
  
  std::cout << ((pm.is_valid()) ? "map valid!" : "map invalid!") << std::endl
            << std::endl;

  Point_2 p1(2, 2);
  if (!pm.is_point_in_face(p1, pm.unbounded_face())) {
    std::cerr << "Error: failed to determine point1 in unbounded face!"
              << std::endl;
    return 1;
  }
  
  if (e[0]->face()->is_unbounded()) {
    std::cerr << "Error: wrong unbounded face!" << std::endl;
    return 1;
  }
  
  Point_2 p2(1, 1);
  if (!pm.is_point_in_face(p2, e[0]->face())) {
    std::cerr << "Error: failed to determine point2 in face!" << std::endl;
    return 1;
  }

  if (pm.is_point_in_face(p1, e[0]->face())) {
    std::cerr << "Error: failed to determine point1 in face!" << std::endl;
    return 1;
  }

  if (pm.is_point_in_face(p2, pm.unbounded_face())) {
    std::cerr << "Error: failed to determine point2 in unbounded face!"
              << std::endl;
    return 1;
  }

  return 0;
}
